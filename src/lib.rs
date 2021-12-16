use std::cmp::Ordering;
use std::fs;

#[derive(Debug)]
pub struct FastaFile {
    pub records: Vec<FastaRecord>,
    pub len: usize,
}

#[derive(Debug)]
pub struct FastaRecord {
    pub id: String,
    pub desc: String,
    pub seq: String,
    pub len: usize,
}

impl FastaFile {
    pub fn new(records: &Vec<FastaRecord>) -> FastaFile {
        FastaFile {
            records: records.to_owned(),
            len: records.len(),
        }
    }

    pub fn single_record_from_file(path: &str) -> FastaRecord {
        let file = fs::read_to_string(path).expect("Couldn't open file");

        let lines: Vec<&str> = file.lines().map(|line| line.trim()).collect();
        let id;
        let desc;
        let seq = lines[1..].join("");

        match lines[0].split_once(&[' ', '|'][..]) {
            Some(id_and_desc) => {
                id = id_and_desc.0.replace(">", "");
                desc = id_and_desc.1;
            }
            None => {
                id = lines[0].replace(">", "");
                desc = "";
            }
        }

        FastaRecord::new(id.as_str(), desc, seq.as_str())
    }

    pub fn parse_multifasta_from_file(path: &str) -> FastaFile {
        let file = fs::read_to_string(path).expect("Couldn't open file");

        let mut records: Vec<FastaRecord> = Vec::new();

        let mut id;
        let mut desc;
        let mut seq = String::new();

        let mut fasta_lines = file.lines().map(|line| line.trim());

        let first_header = fasta_lines.next().expect("Couldn't extract first line");

        match first_header.split_once(&[' ', '|'][..]) {
            Some(id_and_desc) => {
                id = id_and_desc.0.replace(">", "");
                desc = id_and_desc.1;
            }
            None => {
                id = first_header.replace(">", "");
                desc = "";
            }
        }

        for line in fasta_lines {
            if line.is_empty() {
                continue;
            } else if line.as_bytes()[0] == b'>' {
                records.push(FastaRecord::new(id.as_str(), desc, seq.as_str()));
                seq = String::new();

                match line.split_once(&[' ', '|'][..]) {
                    Some(id_and_desc) => {
                        id = id_and_desc.0.replace(">", "");
                        desc = id_and_desc.1;
                    }
                    None => {
                        id = line.replace(">", "");
                        desc = "";
                    }
                }
            } else {
                seq += line;
            }
        }

        // Make sure last record goes in
        records.push(FastaRecord::new(id.as_str(), desc, seq.as_str()));

        FastaFile::new(&records)
    }
}

impl FastaRecord {
    pub fn new(id: &str, desc: &str, seq: &str) -> FastaRecord {
        FastaRecord {
            id: id.to_string(),
            desc: desc.to_string(),
            seq: seq.to_string(),
            len: seq.len(),
        }
    }

    pub fn gc_content(&self) -> f32 {
        let gc_count = self
            .seq
            .chars()
            .filter(|&nuc| nuc == 'C' || nuc == 'G')
            .count() as f32;

        gc_count / self.len as f32
    }
}

impl PartialEq for FastaFile {
    fn eq(&self, other: &Self) -> bool {
        if self.len != other.len {
            return false;
        }

        let mut self_records = self.records.clone();
        let mut other_records = other.records.clone();

        self_records.sort();
        other_records.sort();

        for i in 0..self.len {
            if self_records[i] != other_records[i] {
                return false;
            }
        }

        true
    }
}

impl PartialEq for FastaRecord {
    fn eq(&self, other: &FastaRecord) -> bool {
        self.id == other.id && self.seq == other.seq && self.len == other.len
    }
}

impl Eq for FastaRecord {}

impl Clone for FastaRecord {
    fn clone(&self) -> Self {
        FastaRecord::new(self.id.as_str(), self.desc.as_str(), self.seq.as_str())
    }
}

impl PartialOrd for FastaRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for FastaRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.id.cmp(&other.id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::{self, ApproxEq};

    #[test]
    fn creates_fasta_record() {
        let file = fs::read_to_string("./test.fas").expect("Couldn't read file");
        let expected = FastaRecord {
            id: String::from("test"),
            desc: String::from(""),
            seq: String::from("ACTGTGAC"),
            len: 8,
        };
        let parsed: Vec<&str> = file
            .split_terminator(&['\n', '\r'][..])
            .map(|line| line.trim())
            .collect();

        let id = parsed[0].replace(">", "");
        let seq = parsed[1];

        let test = FastaRecord::new(id.as_str(), "", seq);

        assert_eq!(expected, test);
    }

    #[test]
    fn loads_single_record_from_fasta() {
        let path = "./test.fas";
        let expected = FastaRecord {
            id: String::from("test"),
            desc: String::from(""),
            seq: String::from("ACTGTGAC"),
            len: 8,
        };
        let parsed = FastaFile::single_record_from_file(path);

        assert_eq!(expected, parsed);
    }

    #[test]
    fn loads_multifasta() {
        let path = "./test_multiple.fas";
        let parsed = FastaFile::parse_multifasta_from_file(path);

        let expected = FastaFile::new(&vec![
            FastaRecord::new("test1", "desc1", "ACTG"),
            FastaRecord::new("test2", "desc2", "TGCA"),
            FastaRecord::new("test3", "desc3", "TTGAGA"),
        ]);

        assert_eq!(parsed, expected);
    }

    #[test]
    fn loads_multiline_seq_from_multifasta() {
        let path = "./test_multiline.fas";
        let parsed = FastaFile::parse_multifasta_from_file(path);

        let expected = FastaFile::new(&vec![
            FastaRecord::new("test1", "multiline1", "ACTGTGCATGAC"),
            FastaRecord::new("test2", "multiline2 long description", "ACTGTGACTGTGAC"),
        ]);

        assert_eq!(parsed, expected);
    }

    #[test]

    fn calculates_gc_content() {
        let path = "./test_gc.fas";
        let record = FastaFile::single_record_from_file(path);
        let target = 60.919540;

        assert!(target.approx_eq(record.gc_content() * 100.0, (0.001, 2)));
    }
}
