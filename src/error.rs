use std::fmt;

use gdal_typed_rasterband::typed_rasterband::TypeError;

#[derive(Debug)]
pub struct Error {
    msg: String,
}

impl Error {
    pub fn from_string(msg: String) -> Error {
        Error { msg }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.msg)
    }
}

impl std::error::Error for Error {
    fn description(&self) -> &str {
        &self.msg
    }

    fn cause(&self) -> Option<&std::error::Error> {
        None
    }
}

impl std::convert::From<TypeError> for Error {
    fn from(_error: TypeError) -> Error {
        Error {
            msg: "band TypeError".to_string(),
        }
    }
}

impl std::convert::From<gdal::errors::Error> for Error {
    fn from(error: gdal::errors::Error) -> Error {
        Error {
            msg: format!("GdalError {}", error),
        }
    }
}
