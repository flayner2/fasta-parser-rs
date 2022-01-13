use phf::{phf_map, Map};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fs;

pub const CODON_TABLE: Map<&str, &str> = phf_map! {
    "UUU" => "F",
    "CUU" => "L",
    "AUU" => "I",
    "GUU" => "V",
    "UUC" => "F",
    "CUC" => "L",
    "AUC" => "I",
    "GUC" => "V",
    "UUA" => "L",
    "CUA" => "L",
    "AUA" => "I",
    "GUA" => "V",
    "UUG" => "L",
    "CUG" => "L",
    "AUG" => "M",
    "GUG" => "V",
    "UCU" => "S",
    "CCU" => "P",
    "ACU" => "T",
    "GCU" => "A",
    "UCC" => "S",
    "CCC" => "P",
    "ACC" => "T",
    "GCC" => "A",
    "UCA" => "S",
    "CCA" => "P",
    "ACA" => "T",
    "GCA" => "A",
    "UCG" => "S",
    "CCG" => "P",
    "ACG" => "T",
    "GCG" => "A",
    "UAU" => "Y",
    "CAU" => "H",
    "AAU" => "N",
    "GAU" => "D",
    "UAC" => "Y",
    "CAC" => "H",
    "AAC" => "N",
    "GAC" => "D",
    "UAA" => "*",
    "CAA" => "Q",
    "AAA" => "K",
    "GAA" => "E",
    "UAG" => "*",
    "CAG" => "Q",
    "AAG" => "K",
    "GAG" => "E",
    "UGU" => "C",
    "CGU" => "R",
    "AGU" => "S",
    "GGU" => "G",
    "UGC" => "C",
    "CGC" => "R",
    "AGC" => "S",
    "GGC" => "G",
    "UGA" => "*",
    "CGA" => "R",
    "AGA" => "R",
    "GGA" => "G",
    "UGG" => "W",
    "CGG" => "R",
    "AGG" => "R",
    "GGG" => "G",
};

fn get_id_and_description(header: &str) -> (String, String) {
    let (id, desc) = header
        .split_once(&[' ', '|'][..])
        .unwrap_or_else(|| (header, ""));

    let id = id.replace(">", "");

    (id, desc.to_owned())
}

#[derive(Debug, Eq, Clone)]
pub enum SeqType {
    Nucleotide,
    Protein,
}

#[derive(Debug, Eq, Clone)]
pub struct FastaFile {
    pub records: BTreeSet<FastaRecord>,
}

#[derive(Debug, Eq, Clone)]
pub struct FastaRecord {
    pub id: String,
    pub desc: String,
    pub seq: String,
    pub seq_type: SeqType,
}

impl PartialEq for SeqType {
    fn eq(&self, other: &Self) -> bool {
        self == other
    }
}

impl FastaFile {
    pub fn new(records: BTreeSet<FastaRecord>) -> FastaFile {
        FastaFile {
            records: records.to_owned(),
        }
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn single_record_from_file(
        path: &str,
        nucleotide: Option<bool>,
    ) -> Result<FastaRecord, std::io::Error> {
        let file = fs::read_to_string(path)?;

        let lines: Vec<&str> = file.lines().map(|line| line.trim()).collect();
        let seq = lines[1..].join("");
        let seq_type = FastaRecord::infer_type(&seq, nucleotide);

        let (id, desc) = get_id_and_description(lines[0]);

        Ok(FastaRecord::new(id, desc, seq, seq_type))
    }

    pub fn parse_multifasta_from_file(
        path: &str,
        nucleotide: Option<bool>,
    ) -> Result<FastaFile, std::io::Error> {
        let file = fs::read_to_string(path)?;

        let mut records: BTreeSet<FastaRecord> = BTreeSet::new();

        let mut seq_type;
        let mut seq = String::new();

        let mut fasta_lines = file.lines().map(|line| line.trim());

        let first_header = fasta_lines.next();

        let (mut id, mut desc) = match first_header {
            Some(first_header) => get_id_and_description(first_header),
            None => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Input file is empty, corrupted or inaccessible.",
                ))
            }
        };

        for line in fasta_lines.filter(|line| !line.is_empty()) {
            if line.starts_with('>') {
                seq_type = FastaRecord::infer_type(&seq, nucleotide);
                records.insert(FastaRecord::new(id, desc, seq, seq_type));
                seq = String::new();

                let id_and_desc = get_id_and_description(line);

                id = id_and_desc.0;
                desc = id_and_desc.1;
            } else {
                seq += line;
            }
        }

        // Make sure last record goes in
        seq_type = FastaRecord::infer_type(&seq, nucleotide);
        records.insert(FastaRecord::new(id, desc, seq, seq_type));

        Ok(FastaFile::new(records))
    }
}

impl FastaRecord {
    pub fn new(id: String, desc: String, seq: String, seq_type: SeqType) -> FastaRecord {
        FastaRecord {
            id,
            desc,
            seq,
            seq_type,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn gc_content(&self) -> f32 {
        let gc_count = self
            .seq
            .bytes()
            .filter(|&nuc| nuc == b'C' || nuc == b'G')
            .count() as f32;

        gc_count / self.len() as f32
    }

    pub fn transcribe(&self) -> String {
        self.seq.clone().replace("T", "U")
    }

    pub fn translate(&self, should_transcribe: bool) -> String {
        let seq;

        if should_transcribe {
            seq = self.transcribe();
        } else {
            seq = self.seq.clone();
        }

        match self.seq_type {
            SeqType::Protein => seq,
            SeqType::Nucleotide => {
                let mut translated = String::new();

                for codon in seq.as_bytes().chunks(3) {
                    match CODON_TABLE.get(std::str::from_utf8(codon).unwrap()) {
                        None => continue,
                        Some(aminoacid) => {
                            if *aminoacid == "*" {
                                continue;
                            }

                            translated += aminoacid;
                        }
                    }
                }

                translated
            }
        }
    }

    fn infer_type(seq: &str, nucleotide: Option<bool>) -> SeqType {
        match nucleotide {
            Some(is_nuc) => {
                if is_nuc {
                    return SeqType::Nucleotide;
                } else {
                    return SeqType::Protein;
                }
            }
            None => {
                if seq.to_uppercase().contains(
                    &[
                        'F', 'L', 'I', 'V', 'S', 'P', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'W', 'R',
                    ][..],
                ) {
                    return SeqType::Protein;
                } else {
                    return SeqType::Nucleotide;
                }
            }
        }
    }
}

impl PartialEq for FastaFile {
    fn eq(&self, other: &Self) -> bool {
        self.records == other.records
    }
}

impl PartialEq for FastaRecord {
    fn eq(&self, other: &FastaRecord) -> bool {
        self.id == other.id && self.seq == other.seq && self.len() == other.len()
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
            seq_type: SeqType::Nucleotide,
        };
        let parsed: Vec<&str> = file
            .split_terminator(&['\n', '\r'][..])
            .map(|line| line.trim())
            .collect();

        let id = parsed[0].replace(">", "");
        let seq = parsed[1];

        let test = FastaRecord::new(id, String::new(), seq.to_owned(), SeqType::Nucleotide);

        assert_eq!(expected, test);
    }

    #[test]
    fn loads_single_record_from_fasta() -> Result<(), std::io::Error> {
        let path = "./test.fas";
        let expected = FastaRecord {
            id: String::from("test"),
            desc: String::from(""),
            seq: String::from("ACTGTGAC"),
            seq_type: SeqType::Nucleotide,
        };
        let parsed = FastaFile::single_record_from_file(path, None)?;

        assert_eq!(expected, parsed);

        Ok(())
    }

    #[test]
    fn loads_multifasta() -> Result<(), std::io::Error> {
        let path = "./test_multiple.fas";
        let parsed = FastaFile::parse_multifasta_from_file(path, None)?;

        let expected = FastaFile::new(BTreeSet::from([
            FastaRecord::new(
                "test1".to_owned(),
                "desc1".to_owned(),
                "ACTG".to_owned(),
                SeqType::Nucleotide,
            ),
            FastaRecord::new(
                "test2".to_owned(),
                "desc2".to_owned(),
                "TGCA".to_owned(),
                SeqType::Nucleotide,
            ),
            FastaRecord::new(
                "test3".to_owned(),
                "desc3".to_owned(),
                "TTGAGA".to_owned(),
                SeqType::Nucleotide,
            ),
        ]));

        assert_eq!(parsed, expected);

        Ok(())
    }

    #[test]
    fn loads_multiline_seq_from_multifasta() -> Result<(), std::io::Error> {
        let path = "./test_multiline.fas";
        let parsed = FastaFile::parse_multifasta_from_file(path, None)?;

        let expected = FastaFile::new(BTreeSet::from([
            FastaRecord::new(
                "test1".to_owned(),
                "multiline1".to_owned(),
                "ACTGTGCATGAC".to_owned(),
                SeqType::Nucleotide,
            ),
            FastaRecord::new(
                "test2".to_owned(),
                "multiline2 long description".to_owned(),
                "ACTGTGACTGTGAC".to_owned(),
                SeqType::Nucleotide,
            ),
        ]));

        assert_eq!(parsed, expected);

        Ok(())
    }

    #[test]
    fn calculates_gc_content() -> Result<(), std::io::Error> {
        let path = "./test_gc.fas";
        let record = FastaFile::single_record_from_file(path, None)?;
        let target = 60.919540;

        assert!(target.approx_eq(record.gc_content() * 100.0, (0.001, 2)));

        Ok(())
    }

    #[test]
    fn translates() {
        let rna_seq = FastaRecord::new(
            "rna_seq".to_owned(),
            "".to_owned(),
            "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA".to_owned(),
            SeqType::Nucleotide,
        );
        let expected = "MAMAPRTEINSTRING";
        let translated = rna_seq.translate(false);

        assert_eq!(translated, expected);
    }

    #[test]
    fn translates_long() -> Result<(), std::io::Error> {
        let to_translate =
            FastaFile::single_record_from_file("./test_translate_long_rna.fas", None)?;
        let protein = FastaFile::single_record_from_file("./test_translate_long_prot.fas", None)?;

        assert_eq!(to_translate.translate(false), protein.seq);

        Ok(())
    }
}
