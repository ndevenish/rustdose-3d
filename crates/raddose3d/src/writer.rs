//! Writer convenience types for directing output to different destinations.
//!
//! Port of Java's `WriterConsole`, `WriterFile`, `WriterString`, `WriterMultiple`.
//!
//! All types produce `Box<dyn Write + Send + Sync>` which is what the Output modules accept.

use std::fs::File;
use std::io::{self, BufWriter, Cursor, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Open a buffered file for writing.
///
/// Port of Java `WriterFile`.
pub fn file_writer(path: impl AsRef<Path>) -> io::Result<Box<dyn Write + Send + Sync>> {
    let f = File::create(path)?;
    Ok(Box::new(BufWriter::new(f)))
}

/// Write to stdout.
///
/// Port of Java `WriterConsole`.
pub fn stdout_writer() -> Box<dyn Write + Send + Sync> {
    Box::new(io::stdout())
}

/// Write to stderr.
pub fn stderr_writer() -> Box<dyn Write + Send + Sync> {
    Box::new(io::stderr())
}

/// Accumulate all written bytes into an in-memory buffer.
///
/// Port of Java `WriterString`.
/// Call `.into_string()` after use to retrieve the content.
pub struct StringWriter {
    buf: Cursor<Vec<u8>>,
}

impl StringWriter {
    pub fn new() -> Self {
        StringWriter {
            buf: Cursor::new(Vec::new()),
        }
    }

    /// Consume and return the buffered UTF-8 string.
    pub fn into_string(self) -> String {
        String::from_utf8(self.buf.into_inner()).unwrap_or_default()
    }
}

impl Default for StringWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl Write for StringWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.buf.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Forwards all writes to multiple inner writers simultaneously.
///
/// Port of Java `WriterMultiple`.
///
/// Wraps inner writers in `Arc<Mutex<>>` so they can be shared across threads and
/// recovered after use.
#[derive(Clone)]
pub struct TeeWriter {
    writers: Vec<Arc<Mutex<Box<dyn Write + Send + Sync>>>>,
}

impl TeeWriter {
    /// Create a tee writer from a list of boxed writers.
    pub fn new(writers: Vec<Box<dyn Write + Send + Sync>>) -> Self {
        TeeWriter {
            writers: writers
                .into_iter()
                .map(|w| Arc::new(Mutex::new(w)))
                .collect(),
        }
    }
}

impl Write for TeeWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        for w in &self.writers {
            let _ = w.lock().unwrap().write_all(buf);
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        for w in &self.writers {
            let _ = w.lock().unwrap().flush();
        }
        Ok(())
    }
}

/// Discard all output (null device).
pub struct NullWriter;

impl Write for NullWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Construct a boxed `StringWriter` that can be passed to an Output module.
/// The inner `Arc<Mutex<StringWriter>>` can be cloned to read the result later.
pub fn shared_string_writer() -> (Box<dyn Write + Send + Sync>, Arc<Mutex<Vec<u8>>>) {
    let buf: Arc<Mutex<Vec<u8>>> = Arc::new(Mutex::new(Vec::new()));
    let writer = SharedVecWriter(Arc::clone(&buf));
    (Box::new(writer), buf)
}

/// Internal: a `Write` impl backed by a shared `Vec<u8>`.
struct SharedVecWriter(Arc<Mutex<Vec<u8>>>);

impl Write for SharedVecWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

// SharedVecWriter is Send + Sync because Arc<Mutex<_>> is.
unsafe impl Send for SharedVecWriter {}
unsafe impl Sync for SharedVecWriter {}
