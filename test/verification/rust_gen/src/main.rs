use serde::Serialize;
use std::fs::File;
use std::io::Write;
use rand::Rng;
use sassy::{Searcher, profiles::Iupac, Strand};

// --- The Structures for JSON ---
#[derive(Serialize)]
struct TestCase {
    input: InputState,
    output: OutputState,
}

#[derive(Serialize)]
struct InputState {
    hp_in: u64,
    hm_in: u64,
    vp: u64,
    vm: u64,
    eq: u64,
}

#[derive(Serialize)]
struct OutputState {
    hp_out: u64,
    hm_out: u64,
    vp: u64,
    vm: u64,
}

// --- Search Test Types ---
#[derive(Serialize)]
struct SearchTestCase {
    input: SearchInput,
    output: Vec<SearchMatch>,
}

#[derive(Serialize)]
struct SearchInput {
    pattern: String,
    text: String,
    k: usize,
}

#[derive(Serialize)]
struct SearchMatch {
    text_end: usize,
    cost: usize,
}

// --- The Logic Kernel ---
// This is a SCALAR adaptation of your compute_block_simd from bitpacking.rs.
// Since SIMD operations are just parallel scalar operations, running this 
// on u64 is mathematically identical to running it on one lane of AVX2.
fn compute_block_scalar(hp0: &mut u64, hm0: &mut u64, vp: &mut u64, vm: &mut u64, eq: u64) {
    let vx = eq | *vm;
    
    // In Rust code: let eq = eq | *hm0;
    // We rename to eq_c to avoid shadowing confusion in this snippet
    let eq_c = eq | *hm0; 

    // The Myers Folding Logic
    // Rust: (((eq & *vp) + *vp) ^ *vp) | eq
    // Critical: wrapping_add prevents panic on overflow, simulating CPU generic add behavior
    let term = (eq_c & *vp).wrapping_add(*vp); 
    let hx = (term ^ *vp) | eq_c;

    let hp = *vm | !(hx | *vp);
    let hm = *vp & hx;

    // Carries out (to next block)
    // Rust: let hpw = hp >> (u64::BITS - 1);
    let hpw = hp >> 63;
    let hmw = hm >> 63;

    // Carries in (for next vertical step)
    // Rust: let hp = (hp << 1) | *hp0;
    let hp_next = (hp << 1) | *hp0;
    let hm_next = (hm << 1) | *hm0;

    // Update state
    *hp0 = hpw;
    *hm0 = hmw;
    
    let tmp = vx | hp_next;
    *vp = hm_next | !tmp;
    *vm = hp_next & vx;
}

fn generate_block_tests(rng: &mut impl Rng) {
    let mut tests = Vec::new();

    for _ in 0..10_000 {
        // Generate random inputs
        let hp_in: u64 = rng.gen_range(0..2); 
        let hm_in: u64 = rng.gen_range(0..2);
        
        let vp: u64 = rng.gen();
        let vm: u64 = rng.gen();
        let eq: u64 = rng.gen();

        // Clone for running the function
        let mut hp_run = hp_in;
        let mut hm_run = hm_in;
        let mut vp_run = vp;
        let mut vm_run = vm;

        // Execute
        compute_block_scalar(&mut hp_run, &mut hm_run, &mut vp_run, &mut vm_run, eq);

        tests.push(TestCase {
            input: InputState { hp_in, hm_in, vp, vm, eq },
            output: OutputState {
                hp_out: hp_run,
                hm_out: hm_run,
                vp: vp_run,
                vm: vm_run,
            },
        });
    }

    // Save to file
    let file = File::create("../test_vectors.json").expect("Unable to create file");
    serde_json::to_writer_pretty(file, &tests).expect("Unable to write JSON");
    println!("Generated 10,000 test vectors in test_vectors.json");
}

fn generate_search_tests(rng: &mut impl Rng) {
    let mut searcher = Searcher::<Iupac>::new_fwd();
    let mut tests = Vec::new();

    // 1000 random search cases
    for i in 0..1000 {
        // Random Lengths
        let m = rng.gen_range(10..40);
        let n = rng.gen_range(200..1000); // Check enough length to trigger chunking
        let k = rng.gen_range(0..6);

        // Random DNA
        let pattern: String = (0..m).map(|_| {
            match rng.gen_range(0..4) {
                0 => 'A', 1 => 'C', 2 => 'G', _ => 'T'
            }
        }).collect();
        
        let mut text: String = (0..n).map(|_| {
            match rng.gen_range(0..4) {
                0 => 'A', 1 => 'C', 2 => 'G', _ => 'T'
            }
        }).collect();

        // Inject match 50% of the time, ensure we have some hits
        if i % 2 == 0 {
             let pos = rng.gen_range(0..(n - m));
             // Replace slice in string... logic is a bit annoying in Rust string,
             // convert to Vec<u8> first used below anyway.
        }

        let p_bytes = pattern.as_bytes();
        let mut t_bytes = text.as_bytes().to_vec();

        if i % 2 == 0 {
             let pos = rng.gen_range(0..(n - m));
             t_bytes[pos..pos+m].copy_from_slice(p_bytes);
             // Update text string
             text = String::from_utf8(t_bytes.clone()).unwrap();
        }

        // Run Sassy
        let matches = searcher.search(p_bytes, &t_bytes, k);

        // Convert Matches
        let mut out_matches = Vec::new();
        for m in matches {
            out_matches.push(SearchMatch {
                text_end: m.text_end, // 0-based exclusive end
                cost: m.cost as usize,
            });
        }
        
        // Sort for deterministic validation order
        out_matches.sort_by(|a, b| a.text_end.cmp(&b.text_end));

        tests.push(SearchTestCase {
            input: SearchInput { pattern, text, k },
            output: out_matches,
        });
    }

    let file = File::create("../test_vectors_search.json").expect("Unable to create file");
    serde_json::to_writer_pretty(file, &tests).expect("Unable to write JSON");
    println!("Generated 1,000 search vectors in test_vectors_search.json");
}

fn main() {
    let mut rng = rand::thread_rng();
    generate_block_tests(&mut rng);
    generate_search_tests(&mut rng);
}
