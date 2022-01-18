///Shows simple usage of raw metis routine
///Uses a graph copied directly from METIS manual

use sparse_qr_rust::metis::METIS_ComputeVertexSeparator;

fn main() {
    let xadj : Vec<i64> = vec![0,2,5,8,11,13,16,20,24,28,31,33,36,39,42,44];
    let adjncy : Vec<i64> = vec![1,5,0,2,6,1,3,7,2,4,8,3,9,0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14,5,11,6,10,12,7,11,13,8,12,14,9,13];
    let nvtxs : i64 = (xadj.len()-1) as i64;
    let mut part : Vec<i64> = vec![0;nvtxs as usize];
    let mut sepsize : i64 = 0;

    unsafe{
        METIS_ComputeVertexSeparator(&nvtxs,xadj.as_ptr(),adjncy.as_ptr(),std::ptr::null(),std::ptr::null(),&mut sepsize,part.as_mut_ptr());
    }


    println!("sepsize = {:?}",sepsize);
    println!("part =  {:?}",part);

}
