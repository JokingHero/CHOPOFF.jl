use serde::Serialize;
use std::fs::File;
use std::path::Path;

#[derive(Serialize)]
struct EncodingTestCase {
    text: String,
    encoded: Vec<u64>, 
    bases: Vec<u8>, 
}

// --- IUPAC CONSTANTS (Copied from sassy/profiles/iupac.rs) ---

#[rustfmt::skip]
const IUPAC_CODE: [u8; 32] = {
    let mut t = [255u8; 32];
    const A: u8 = 1 << 0;
    const C: u8 = 1 << 1;
    const T: u8 = 1 << 2;
    const G: u8 = 1 << 3;

    t[b'A' as usize & 0x1F] = A;
    t[b'C' as usize & 0x1F] = C;
    t[b'T' as usize & 0x1F] = T;
    t[b'U' as usize & 0x1F] = T;
    t[b'G' as usize & 0x1F] = G;
    t[b'N' as usize & 0x1F] = A|C|T|G;
    
    t[b'R' as usize & 0x1F] = A|G;
    t[b'Y' as usize & 0x1F] = C|T;
    t[b'S' as usize & 0x1F] = G|C;
    t[b'W' as usize & 0x1F] = A|T;
    t[b'K' as usize & 0x1F] = G|T;
    t[b'M' as usize & 0x1F] = A|C;
    t[b'B' as usize & 0x1F] = C|G|T;
    t[b'D' as usize & 0x1F] = A|G|T;
    t[b'H' as usize & 0x1F] = A|C|T;
    t[b'V' as usize & 0x1F] = A|C|G;
    
    t[b'X' as usize & 0x1F] = 0;
    
    t
};

fn get_encoded(c: u8) -> u8 {
    IUPAC_CODE[(c & 0x1F) as usize]
}

// --- LOGIC REPLICATION ---

// This mimics `encode_ref` but using scalar loop instead of SIMD shuffles.
// Since the SIMD shuffle logic just implements this lookup table in parallel,
// this scalar version is mathematically equivalent for verification.
fn encode_ref_scalar(text: &[u8; 64], bases: &[u8]) -> Vec<u64> {
    let mut masks = vec![0u64; bases.len()];

    for (j, &base) in bases.iter().enumerate() {
        let base_mask = get_encoded(base);
        let mut result = 0u64;

        // Iterate 0..64. 
        for i in 0..64 {
            let text_char = text[i];
            let text_mask = get_encoded(text_char);
            
            // Check match
            if (base_mask & text_mask) > 0 {
                result |= 1u64 << i;
            }
        }
        masks[j] = result;
    }
    masks
}

fn main() {
    // 1. Define sequence
    let seq_str = "ACGTacgtNnRrYy"; 
    let seq_bytes = seq_str.as_bytes();
    
    // 2. Pad to 64 bytes
    let mut block = [b'X'; 64];
    for (i, &b) in seq_bytes.iter().enumerate() {
        block[i] = b;
    }

    // 3. Define Bases we want to encode (Standard DNA)
    let bases = vec![b'A', b'C', b'T', b'G'];

    // 4. Run Encoding
    let encoded = encode_ref_scalar(&block, &bases);

    // 5. Output
    let test_case = EncodingTestCase {
        text: seq_str.to_string(),
        encoded,
        bases,
    };

    let path = Path::new("../encoding_test.json");
    let file = File::create(&path).expect("Could not create file");
    serde_json::to_writer_pretty(file, &test_case).expect("Could not write JSON");
    
    println!("SUCCESS: Generated {:?} with sequence: {}", path, seq_str);
}
