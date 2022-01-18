use sparse_qr_rust::dtree::DissectionTree;
use sparse_qr_rust::numeric::SparseQR;
use sparse_qr_rust::gallery::laplace3d;



fn main(){
    use std::env;
    let args: Vec<String> = env::args().collect();
    if args.len()<3{
        panic!("usage: ./laplace3d_factorize extent maxnodes");
    }
    let n=args[1].parse::<usize>().unwrap();
    let maxnodes=args[2].parse::<usize>().unwrap();

    let mx = n;
    let my = n;
    let mz = n;
    let m = mx*my*mz;
    type F=f64;
    //Assemble matrix
    let lap=laplace3d::<F>(mx,my,mz);
    //Get squared graph
    let g = lap.to_metis_graph();
    let g2 = g.square();
    //Symbolic factorization
    let dtree = DissectionTree::new(&g2,maxnodes);
    //Numeric factorization, consumes input tree
    let mut fact = SparseQR::<F>::new(dtree,&lap);

    //Make some data to solve with
    let soln : Vec<F> = (0..m).map(|x| (x as f64).sin()+2.0).collect();
    let mut b = soln.clone();
    lap.eval(&soln,&mut b);

    //Now go back and solve
    fact.solve(&mut b);

    //Result should be close to `soln`
    let mut errs : Vec<F> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).abs()/x.abs()).collect();
    errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
    let maxerr = errs.last().unwrap();
    println!("Maximum relative error: {}",maxerr);
    println!("Bytes by level");
    fact.print_level_bytes();
    println!("Timings by level");
    fact.print_level_times();
}
