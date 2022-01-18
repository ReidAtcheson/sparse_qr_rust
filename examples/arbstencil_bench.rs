use sparse_qr_rust::dtree::DissectionTree;
use sparse_qr_rust::numeric::SparseQR;
use sparse_qr_rust::gallery::{arbstencil2d,arbstencil3d};



fn main(){
    use std::env;
    let args: Vec<String> = env::args().collect();
    if args.len()<5{
        panic!("usage: ./laplace3d_factorize dtype dim extent maxnodes");
    }
    let dtype=args[1].clone();
    let dim=args[2].parse::<usize>().unwrap();
    let n=args[3].parse::<usize>().unwrap();
    let maxnodes=args[4].parse::<usize>().unwrap();

    assert!(dtype=="f32" || dtype=="f64" || dtype=="c32" || dtype=="c64");
    assert!(dim==2 || dim==3);
    assert!(n>0);
    assert!(maxnodes>0);
    if dtype=="f64"{
        type F=f64;
        let (mat,m) = {
            let mx = n;
            let my = n;
            let mz = n;
            let m = if dim==2 { mx*my } else { mx*my*mz };
            let mat = if dim==2 {arbstencil2d::<F>(mx,my,&[-5.0,1.1,0.9,0.8,1.5])} else {arbstencil3d::<F>(mx,my,mz,&[-8.0,1.1,0.9,1.5,1.7,0.9,0.5])};
            (mat,m)
        };
        //Get squared graph
        let g = mat.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&mat);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).map(|x| (x as F).sin()+2.0).collect();
        let mut b = vec![0.0;soln.len()];
        mat.eval(&soln,&mut b);

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

    if dtype=="f32"{
        type F=f32;
        let (mat,m) = {
            let mx = n;
            let my = n;
            let mz = n;
            let m = if dim==2 { mx*my } else { mx*my*mz };
            let mat = if dim==2 {arbstencil2d::<F>(mx,my,&[-5.0,1.1,0.9,0.8,1.5])} else {arbstencil3d::<F>(mx,my,mz,&[-8.0,1.1,0.9,1.5,1.7,0.9,0.5])};
            (mat,m)
        };
        //Get squared graph
        let g = mat.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&mat);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).map(|x| (x as F).sin()+2.0).collect();
        let mut b = vec![0.0;soln.len()];
        mat.eval(&soln,&mut b);

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

    if dtype=="c32"{
        use num_complex::Complex;
        type R=f32;
        type F=Complex<R>;
        let (mat,m) = {
            let mx = n;
            let my = n;
            let mz = n;
            let m = if dim==2 { mx*my } else { mx*my*mz };
            let mat = if dim==2 {
                arbstencil2d::<F>(mx,my,&[F::new(-8.0,3.0),F::new(-1.1,0.3),F::new(0.9,1.2),F::new(1.2,1.1),F::new(0.4,0.3)])
            } 
            else {
                arbstencil3d::<F>(mx,my,mz,&[F::new(-8.0,3.0),F::new(-1.1,0.3),F::new(0.9,1.2),F::new(1.2,1.1),F::new(0.4,0.3),F::new(1.2,-1.1),F::new(0.2,0.5)])
            };
            (mat,m)
        };
        //Get squared graph
        let g = mat.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&mat);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).zip(0..m).map(|(xr,xi)| {
            F::new((xr as R).sin(),(xi as R).cos())
        }).collect();


        let mut b = vec![F::new(0.0,0.0);soln.len()];
        mat.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<R> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).norm()/x.norm()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        println!("Maximum relative error: {}",maxerr);
        println!("Bytes by level");
        fact.print_level_bytes();
        println!("Timings by level");
        fact.print_level_times();
    }

    if dtype=="c64"{
        use num_complex::Complex;
        type R=f64;
        type F=Complex<R>;
        let (mat,m) = {
            let mx = n;
            let my = n;
            let mz = n;
            let m = if dim==2 { mx*my } else { mx*my*mz };
            let mat = if dim==2 {
                arbstencil2d::<F>(mx,my,&[F::new(-8.0,3.0),F::new(-1.1,0.3),F::new(0.9,1.2),F::new(1.2,1.1),F::new(0.4,0.3)])
            } 
            else {
                arbstencil3d::<F>(mx,my,mz,&[F::new(-8.0,3.0),F::new(-1.1,0.3),F::new(0.9,1.2),F::new(1.2,1.1),F::new(0.4,0.3),F::new(1.2,-1.1),F::new(0.2,0.5)])
            };
            (mat,m)
        };
        //Get squared graph
        let g = mat.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&mat);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).zip(0..m).map(|(xr,xi)| {
            F::new((xr as R).sin(),(xi as R).cos())
        }).collect();


        let mut b = vec![F::new(0.0,0.0);soln.len()];
        mat.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<R> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).norm()/x.norm()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        println!("Maximum relative error: {}",maxerr);
        println!("Bytes by level");
        fact.print_level_bytes();
        println!("Timings by level");
        fact.print_level_times();
    }



}
