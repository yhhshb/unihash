//! Minimal streaming FASTA/FASTQ reader (format auto-detected per file,
//! gzip auto-detected by magic bytes when the `gz` feature is enabled).
//! Only sequences are returned; headers and qualities are skipped.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Format {
    Unknown,
    Fasta,
    Fastq,
}

pub struct FastxReader<R: BufRead> {
    inner: R,
    format: Format,
    line: Vec<u8>,
    /// FASTA: a header line for the next record has already been consumed.
    header_pending: bool,
    /// FASTA: sequence bytes are coalesced here (line breaks removed) and
    /// handed out in chunks of `chunk_size` bytes.
    chunk: Vec<u8>,
    chunk_size: usize,
}

/// Default size of the sequence chunks delivered for FASTA records.
pub const DEFAULT_CHUNK_SIZE: usize = 64 << 10;

impl<R: BufRead> FastxReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            format: Format::Unknown,
            line: Vec::new(),
            header_pending: false,
            chunk: Vec::new(),
            chunk_size: DEFAULT_CHUNK_SIZE,
        }
    }

    /// Size of the FASTA sequence chunks passed to
    /// [`next_record_chunks`](Self::next_record_chunks) (default 64 KiB).
    pub fn with_chunk_size(mut self, n: usize) -> Self {
        assert!(n > 0);
        self.chunk_size = n;
        self
    }

    fn read_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        let n = self.inner.read_until(b'\n', &mut self.line)?;
        while matches!(self.line.last(), Some(b'\n' | b'\r')) {
            self.line.pop();
        }
        Ok(n > 0)
    }

    /// Read the next record's sequence into `seq` (cleared first).
    /// Returns `Ok(false)` at end of input. Buffers the whole record; prefer
    /// [`next_record_chunks`](Self::next_record_chunks) for chromosomes.
    pub fn next_record(&mut self, seq: &mut Vec<u8>) -> io::Result<bool> {
        seq.clear();
        self.next_record_chunks(|c| seq.extend_from_slice(c))
    }

    /// Stream the next record's sequence to `f` in chunks, without buffering
    /// the record. FASTA sequence is delivered in chunks of `chunk_size`
    /// bytes (the last one shorter), with line breaks and `\r` removed, so
    /// memory stays bounded even for single-line chromosomes and the
    /// per-chunk cost is paid once per 64 KiB rather than once per line.
    /// FASTQ sequence is delivered as one chunk per record.
    /// Returns `Ok(false)` at end of input.
    pub fn next_record_chunks<F: FnMut(&[u8])>(&mut self, mut f: F) -> io::Result<bool> {
        if self.format == Format::Unknown {
            loop {
                if !self.read_line()? {
                    return Ok(false);
                }
                match self.line.first() {
                    None => continue,
                    Some(b'>') => {
                        self.format = Format::Fasta;
                        self.header_pending = true;
                        break;
                    }
                    Some(b'@') => {
                        self.format = Format::Fastq;
                        // header consumed; read sequence below
                        return self.fastq_body(&mut f);
                    }
                    Some(&c) => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("not FASTA/FASTQ: first character {:?}", c as char),
                        ))
                    }
                }
            }
        }
        match self.format {
            Format::Fasta => self.fasta_record(&mut f),
            Format::Fastq => loop {
                if !self.read_line()? {
                    return Ok(false);
                }
                if self.line.is_empty() {
                    continue;
                }
                if self.line[0] != b'@' {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "malformed FASTQ: expected '@' header",
                    ));
                }
                return self.fastq_body(&mut f);
            },
            Format::Unknown => unreachable!(),
        }
    }

    fn fastq_body<F: FnMut(&[u8])>(&mut self, f: &mut F) -> io::Result<bool> {
        if !self.read_line()? {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ record",
            ));
        }
        f(&self.line);
        // '+' line and quality line
        if !self.read_line()? || !self.read_line()? {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ record",
            ));
        }
        Ok(true)
    }

    fn fasta_record<F: FnMut(&[u8])>(&mut self, f: &mut F) -> io::Result<bool> {
        if !self.header_pending {
            // find next header
            loop {
                if !self.read_line()? {
                    return Ok(false);
                }
                if self.line.first() == Some(&b'>') {
                    break;
                }
            }
        }
        self.header_pending = false;
        self.chunk.clear();
        let mut at_line_start = true;
        loop {
            let buf = self.inner.fill_buf()?;
            if buf.is_empty() {
                break; // EOF ends the record
            }
            if at_line_start && buf[0] == b'>' {
                // Next header: consume it and stop.
                self.read_line()?;
                self.header_pending = true;
                break;
            }
            // Consume up to and including the next newline (or the buffer).
            let (end, newline) = match buf.iter().position(|&c| c == b'\n') {
                Some(p) => (p + 1, true),
                None => (buf.len(), false),
            };
            let body = if newline {
                &buf[..end - 1]
            } else {
                &buf[..end]
            };
            for part in body.split(|&c| c == b'\r') {
                let mut part = part;
                while !part.is_empty() {
                    let room = self.chunk_size - self.chunk.len();
                    let take = room.min(part.len());
                    self.chunk.extend_from_slice(&part[..take]);
                    part = &part[take..];
                    if self.chunk.len() == self.chunk_size {
                        f(&self.chunk);
                        self.chunk.clear();
                    }
                }
            }
            at_line_start = newline;
            self.inner.consume(end);
        }
        if !self.chunk.is_empty() {
            f(&self.chunk);
            self.chunk.clear();
        }
        Ok(true)
    }
}

/// Open a path (or `-` for stdin), transparently decompressing gzip.
pub fn open(path: &Path) -> io::Result<FastxReader<Box<dyn BufRead>>> {
    let raw: Box<dyn Read> = if path.as_os_str() == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| {
            io::Error::new(e.kind(), format!("cannot open {}: {e}", path.display()))
        })?)
    };
    let mut br = BufReader::with_capacity(1 << 20, raw);
    let magic = br.fill_buf()?;
    let is_gz = magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b;
    let reader: Box<dyn BufRead> = if is_gz {
        #[cfg(feature = "gz")]
        {
            Box::new(BufReader::with_capacity(
                1 << 20,
                flate2::read::MultiGzDecoder::new(br),
            ))
        }
        #[cfg(not(feature = "gz"))]
        {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "gzip input requires the `gz` feature",
            ));
        }
    } else {
        Box::new(br)
    };
    Ok(FastxReader::new(reader))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(s: &str) -> Vec<String> {
        let mut r = FastxReader::new(s.as_bytes());
        let mut v = Vec::new();
        let mut buf = Vec::new();
        while r.next_record(&mut buf).unwrap() {
            v.push(String::from_utf8(buf.clone()).unwrap());
        }
        v
    }

    #[test]
    fn fasta() {
        assert_eq!(all(">a\nACG\nTT\n>b\n\nGG\n"), vec!["ACGTT", "GG"]);
        assert_eq!(all(">a\r\nAC\r\n"), vec!["AC"]);
        assert_eq!(all(""), Vec::<String>::new());
    }

    #[test]
    fn chunked_matches_buffered_with_tiny_buffers() {
        let text = ">c1 desc\r\nACGTN\r\nacgt\r\n\r\n>c2\nTTTTTTTTTTTTTTTTTTTT\n>c3\n>c4\nGG";
        for cap in [1usize, 2, 3, 5, 7, 64] {
            let chunk_size = cap.max(2) - 1;
            let mut r = FastxReader::new(io::BufReader::with_capacity(cap, text.as_bytes()))
                .with_chunk_size(chunk_size);
            let mut got = Vec::new();
            loop {
                let mut rec = Vec::new();
                let mut lens = Vec::new();
                if !r
                    .next_record_chunks(|c| {
                        lens.push(c.len());
                        assert!(!c.is_empty() && c.len() <= chunk_size);
                        assert!(!c.contains(&b'\n') && !c.contains(&b'\r'));
                        rec.extend_from_slice(c);
                    })
                    .unwrap()
                {
                    break;
                }
                // Every chunk but the last of a record is full.
                if let Some((_, full)) = lens.split_last() {
                    assert!(full.iter().all(|&n| n == chunk_size), "cap={cap} {lens:?}");
                }
                got.push(String::from_utf8(rec).unwrap());
            }
            assert_eq!(
                got,
                vec!["ACGTNacgt", "TTTTTTTTTTTTTTTTTTTT", "", "GG"],
                "cap={cap}"
            );
        }
    }

    #[test]
    fn fastq() {
        assert_eq!(
            all("@r1\nACGT\n+\nIIII\n@r2\nGG\n+r2\n@@\n"),
            vec!["ACGT", "GG"]
        );
    }
}
