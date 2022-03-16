use crate::dtree::DissectionTree;
use crate::sparse::CSCSparse;
use crate::utility::Lapack;
use std::collections::BTreeMap;
use num_traits::Num;
use num_traits::cast::ToPrimitive;


fn print_diag<F : std::fmt::Display>(nrows : usize, ncols : usize,a : &Vec<F>)->(){
    for i in 0..std::cmp::min(nrows,ncols){
        let j=i;
        print!("{}\n",a[i+nrows*j]);
    }
}

fn gather_matrices<F : Num+Clone>(ncols : usize,data : &BTreeMap<usize,Vec<F>>,ids : &Vec<usize>) -> Vec<F>{
    for k in ids.iter(){
        assert_eq!(data.get(k).unwrap().len() % ncols,0);
    }
    let nrows=ids.iter().map(|k|data.get(k).unwrap().len()/ncols).fold(0,|acc,x|acc+x);
    let mut out = vec![F::zero();nrows*ncols];

    let mut offs=0;
    for k in ids.iter(){
        let m=data.get(k).unwrap();
        let lnrows=m.len()/ncols;
        for (xs,ys) in out.chunks_exact_mut(nrows).zip(m.chunks_exact(lnrows)){
            for (x,y) in (&mut xs[offs..offs+lnrows]).into_iter().zip(ys.iter()){
                *x=y.clone();
            }
        }
        offs+=lnrows;
    }
    out
}
fn scatter_matrices<F : Num+Clone>(ncols : usize,data : &mut BTreeMap<usize,Vec<F>>,ids : &Vec<usize>,stacked : &Vec<F>) -> (){
    for k in ids.iter(){
        assert_eq!(data.get(k).unwrap().len() % ncols,0);
    }
    let nrows=ids.iter().map(|k|data.get(k).unwrap().len()/ncols).fold(0,|acc,x|acc+x);
    assert_eq!(stacked.len(),nrows*ncols);

    let mut offs=0;
    for k in ids.iter(){
        let m=data.get_mut(k).unwrap();
        let lnrows=m.len()/ncols;
        for (xs,ys) in stacked.chunks_exact(nrows).zip(m.chunks_exact_mut(lnrows)){
            for (x,y) in xs[offs..offs+lnrows].into_iter().zip(ys.iter_mut()){
                *y=x.clone();
            }
        }
        offs+=lnrows;
    }
}

fn permute_out<F : Num+Clone>(nrows : usize,b : &[F],p : &[usize]) -> Vec<F>{
    if p.len()>0{
        assert_eq!(b.len() % nrows,0);
        let ncols=b.len()/nrows;
        let out_nrows=p.len();
        let mut out = vec![F::zero();out_nrows*ncols];
        for (xs,ys) in out.chunks_exact_mut(out_nrows).zip(b.chunks_exact(nrows)){
            for (x,pi) in xs.iter_mut().zip(p.iter()){
                *x = ys[*pi].clone();
            }
        }
        out
    }
    else{
        Vec::<F>::new()
    }
}

fn permute_in<F : Num+Clone>(nrows : usize,b : &mut [F],p : &[usize],pb : &[F]) -> (){
    if p.len()>0{
        assert_eq!(b.len() % nrows,0);
        let ncols=b.len()/nrows;
        let out_nrows=p.len();
        assert_eq!(pb.len(),p.len()*ncols);
        for (xs,ys) in pb.chunks_exact(out_nrows).zip(b.chunks_exact_mut(nrows)){
            for (x,pi) in xs.iter().zip(p.iter()){
                ys[*pi]=x.clone();
            }
        }
    }
}






pub struct SparseQR<F>{
    //Tree data
    parents  : Vec<Option<usize>>,
    children : Vec<Option<(usize,usize)>>,
    levels : Vec<Vec<usize>>,
    nodes : Vec<Vec<usize>>,
    //Factorization data
    nrows : usize,
    ncols : usize,

    //Upper triangular blocks for a block of columns.
    //Blocks of columns corresponding to leafs in
    //the nested dissection tree have no upper
    //triangular blocks and so may have zero
    //length
    triu : Vec<Vec<usize>>,
    //Lower triangular blocks for a block of columns.
    //Including diagonal block.
    tril : Vec<Vec<usize>>,

    //Packed storage so that solve phase can proceed
    //with only two calls to LAPACK.
    //Upper triangular block. 
    triu_num : Vec<Vec<F>>,
    //Lower triangular block.
    tril_num : Vec<Vec<F>>,
    //Householder reflector scalar factors
    taus : Vec<Vec<F>>,
    //Timings for each node
    times : Vec<f64>
} 


impl <F : Lapack<F=F>+Num+ToPrimitive+Copy + std::fmt::Display> SparseQR<F>{
    fn panic_if_invalid(&self) -> () {
        let nnodes=self.nodes.len();
        for n in 0..nnodes{
            //Make sure all descendents of node `n`
            //are represented in node `n`s lower triangular
            //part
            let mut stack = vec![n];
            while let Some(node)=stack.pop(){
                let triu=&self.triu[node];
                if let Some((c1,c2)) = self.children[node]{
                    assert!(triu.contains(&c1));
                    assert!(triu.contains(&c2));
                    stack.push(c1);
                    stack.push(c2);
                }
            }

            //Make sure all parents of node `n`
            //are represented in node `n`s upper triangular
            //part
            let tril=&self.tril[n];
            let mut np=n;
            while let Some(p) = self.parents[np]{
                assert!(tril.contains(&p));
                np=p;
            }
        }

        for n in 0..nnodes{
            //Make sure the first element of
            //the `tril` array for a node
            //is equal to that node
            //(so that the diagonal block
            //is correctly presented to LAPACK)
            let tril=&self.tril[n];
            assert_eq!(tril[0],n);
        }
    }
    pub fn new(dtree : DissectionTree, mat : &CSCSparse<F>) -> Self{
        use std::time::{Instant};
        let nnodes=dtree.nodes.len();
        let mut times = vec![0.0;nnodes];
        //Propagate lower triangular parts down tree
        let mut tril : Vec::<Vec<usize>> = vec![Vec::<usize>::new();nnodes];
        for node in 0..nnodes{
            //Lower triangular part contains the diagonal block as well.
            tril[node].push(node);
            //Now propagate down tree
            let mut stack : Vec<usize> = vec![node];
            while let Some(p) = stack.pop(){
                if let Some((c1,c2)) = dtree.children[p]{
                    tril[c1].push(node);
                    tril[c2].push(node);
                    stack.push(c1);
                    stack.push(c2);
                }
            }
        }
        //Propagate upper triangular parts 
        let mut triu : Vec::<Vec<usize>> = vec![Vec::<usize>::new();nnodes];
        for node in 0..nnodes{
            let mut maybep = dtree.parents[node];
            //Loop through the path to the top node
            while let Some(p) = maybep{
                triu[p].push(node);
                maybep=dtree.parents[p];
            }
        }
        //Sort the upper and lower triangular block lists
        //This step isn't strictly necessary but it does
        //make some logic cleaner later.
        //
        //If we sort these lists then we get the following
        //properties:
        //
        //tril[0] is always a diagonal block
        //tril[i+1] is always a block corresponding to one level down from tril[i]
        //triu[i-1] is always a block corresponding to _at most_ one level up from triu[i]
        for tri in tril.iter_mut(){
            tri.sort();
            tri.reverse();
        }
        for tri in triu.iter_mut(){
            tri.sort();
            tri.reverse();
        }

        //Assemble initial matrices into the 
        //temporary matrix list
        let mut temp : Vec<BTreeMap<usize,Vec<F>>> = vec![BTreeMap::<usize,Vec<F>>::new();nnodes];
        for node in 0..nnodes{
            let cols  = &dtree.nodes[node];
            //Assemble lower triangular submatrices
            for tri in tril[node].iter(){
                let rows = &dtree.nodes[*tri];
                let a = mat.slice_copy(rows,cols);
                temp[node].insert(*tri,a);
            }
            //Assemble upper triangular submatrices
            for tri in triu[node].iter(){
                let rows = &dtree.nodes[*tri];
                let a = mat.slice_copy(rows,cols);
                temp[node].insert(*tri,a);
            }
        }
        //Begin factorization.


        //Loop from bottom level of dissection tree to top level
        let mut tril_num : Vec::<Vec<F>> = vec![Vec::<F>::new();nnodes];
        let mut taus : Vec::<Vec<F>> = vec![Vec::<F>::new();nnodes];
        for level in dtree.levels.iter().rev(){
            for node in level.iter().cloned(){
                let start=Instant::now();
                //Gather lower triangular blocks into a stacked matrix
                let ncols=dtree.nodes[node].len();
                let nrows=tril[node].iter().map(|&k|dtree.nodes[k].len()).fold(0,|acc,x|acc+x);
                let nmap=&temp[node];
                tril_num[node]=gather_matrices(ncols,nmap,&tril[node]);


                //First call queries optimal workspace size
                let trilmat=&mut tril_num[node];
                let tau = {
                    let mut tau : Vec<F> = vec![F::zero();ncols];
                    let lwork : i32 = {
                        let mut lwork=-1 as i32;
                        let mut info=0 as i32;
                        let mut work : Vec<F> = vec![F::zero();1];
                        F::xgeqrf(nrows as i32,ncols as i32,trilmat.as_mut_slice(),nrows as i32,tau.as_mut_slice(),work.as_mut_slice(),lwork,&mut info);
                        assert_eq!(info,0);
                        lwork = work[0].to_f64().unwrap().to_i32().unwrap();
                        lwork
                    };
                    //Second call actually performs factorization
                    {
                        let mut info=0 as i32;
                        let mut work  = vec![F::zero();lwork as usize];
                        F::xgeqrf(nrows as i32,ncols as i32,trilmat.as_mut_slice(),nrows as i32,tau.as_mut_slice(),work.as_mut_slice(),lwork,&mut info);
                        /*
                        if check_zero_diagonal(nrows,ncols,&trilmat){
                            println!("{:?}",trilmat.len());
                            println!("{:?}",*level);
                            println!("{:?}",node);
                            println!("(nrows,nocls): ({:?},{:?})",nrows,ncols);
                            print_mat(nrows,ncols,&trilmat);
                            panic!("Zero found on diagonal of an R matrix (input sparse matrix was singular)");
                        }
                        */
                        assert_eq!(info,0);
                    }
                    tau
                };

                //Iterate up parent paths. 
                //Apply the newly computed Q^T to upper triangular blocks
                //for each parent in this path.
                let mut mp=dtree.parents[node];
                while let Some(p) = mp{
                    //Stack upper triangular parts
                    //into packed matrix storage
                    let triu_nrows=nrows;
                    let triu_ncols=dtree.nodes[p].len();
                    let mut triumat = gather_matrices(triu_ncols,&temp[p],&tril[node]);


                    //First call to query optimal workspace size
                    let lwork : i32 = {
                        let side=b'L';
                        let trans=b'T';
                        let mut work = vec![F::zero();1];
                        let mut lwork = -1 as i32;
                        let mut info = 0;
                        F::xmqr(side,trans,triu_nrows as i32,triu_ncols as i32,tau.len() as i32,trilmat.as_slice(),nrows as i32,tau.as_slice(),
                        triumat.as_mut_slice(),triu_nrows as i32,work.as_mut_slice(),lwork,&mut info);
                        assert_eq!(info,0);
                        lwork = work[0].to_f64().unwrap().to_i32().unwrap();
                        lwork
                    };

                    //Now apply Q^T
                    {
                        let side=b'L';
                        let trans=b'T';
                        let mut work = vec![F::zero();lwork as usize];
                        let mut info = 0;
                        assert!(F::no_nans(&trilmat));
                        F::xmqr(side,trans,triu_nrows as i32,triu_ncols as i32,tau.len() as i32,trilmat.as_slice(),nrows as i32,tau.as_slice(),
                        triumat.as_mut_slice(),triu_nrows as i32,work.as_mut_slice(),lwork,&mut info);
                        assert!(F::no_nans(&trilmat));
                        assert_eq!(info,0);
                    };

                    //Unpack matrices
                    scatter_matrices(triu_ncols,&mut temp[p],&tril[node],&triumat);
                    mp=dtree.parents[p];
                }

                //Save tau scalars for the houesholder reflectors
                taus[node]=tau;
                times[node]=start.elapsed().as_secs_f64();
            }
        }

        //Factorization complete but now we need to re-pack the upper
        //triangular matrices
        let mut triu_num  = vec![Vec::<F>::new();nnodes];
        for (k,tri) in triu.iter().enumerate(){
            let ncols=dtree.nodes[k].len();
            triu_num[k]=gather_matrices(ncols,&temp[k],tri);
        }


        let out = SparseQR { parents : dtree.parents, children : dtree.children, levels : dtree.levels, nodes : dtree.nodes,
        nrows : mat.get_nrows(),  ncols : mat.get_ncols(), triu : triu, tril : tril, triu_num : triu_num, tril_num : tril_num,taus : taus,
        times : times
        };

        out.panic_if_invalid();
        out
    }

    pub fn solve(&mut self,inout : &mut [F]) -> (){
        assert!(F::no_nans(inout));
        assert!(inout.len()>0);
        assert!(inout.len() % self.nrows == 0);
        let nrhs = inout.len() / self.nrows;
        //Loop through levels in reverse order
        //and apply each block's Q^T
        for level in self.levels.iter().rev(){
            for node in level.iter().cloned(){
                //Build permutation vector corresponding to lower triangular part
                let perm : Vec<usize> = self.tril[node].iter().map(|&n|self.nodes[n].clone()).flatten().collect();
                let nrows=perm.len();
                //let ncols=self.nodes[node].len();
                //Permute active data into temporary array
                let mut tmp = permute_out(self.nrows,inout,&perm);
                //Apply Q^T to the temporary array
                let qt=&self.tril_num[node];
                let tau=&self.taus[node];
                //First call to query optimal workspace size
                let lwork : i32 = {
                    let side=b'L';
                    let trans=b'T';
                    let mut work = vec![F::zero();1];
                    let mut lwork = -1 as i32;
                    let mut info = 0;
                    F::xmqr(side,trans,nrows as i32,nrhs as i32,tau.len() as i32,qt.as_slice(),nrows as i32,tau.as_slice(),
                    tmp.as_mut_slice(),nrows as i32,work.as_mut_slice(),lwork,&mut info);
                    assert_eq!(info,0);
                    lwork = work[0].to_f64().unwrap().to_i32().unwrap();
                    lwork
                };

                //Now actually apply Q^T
                {
                    let side=b'L';
                    let trans=b'T';
                    let mut work = vec![F::zero();lwork as usize];
                    let mut info = 0;
                    F::xmqr(side,trans,nrows as i32,nrhs as i32,tau.len() as i32,qt.as_slice(),nrows as i32,tau.as_slice(),
                    tmp.as_mut_slice(),nrows as i32,work.as_mut_slice(),lwork,&mut info);
                    assert_eq!(info,0);
                };
                //Permute data from temporary array back into in/out array
                permute_in(self.nrows,inout,&perm,&tmp);
            }
        }
        //Now `inout` has been overwritten with ` Q^T * inout`. We now do back-substitution
        //which is similar to previous loop but in reversed order starting
        //with the top-most node of the nested dissection tree
        for level in self.levels.iter(){
            for node in level.iter().cloned(){
                //Get permutation vector corresponding to upper triangular part
                let perm : Vec<usize> = self.triu[node].iter().map(|&n|self.nodes[n].clone()).flatten().collect();
                //Permutation corresponding to diagonal block
                let permd : Vec<usize> = self.nodes[node].clone();
                let nrows=perm.len();
                let nrowsd=permd.len();
                let ncols=self.nodes[node].len();
                //Permute upper triangular active data into temporary array
                let mut tmp = permute_out(self.nrows,inout,&perm);
                //Permute diagonal block active data into temporary array
                let mut tmpd = permute_out(self.nrows,inout,&permd);


                //Back-substitution of diagonal block
                {
                    let rmat = &self.tril_num[node];
                    let uplo : u8 =b'U';
                    let trans : u8 =b'N';
                    let diag : u8 =b'N';
                    let n=ncols as i32;
                    assert_eq!(rmat.len() % ncols, 0);
                    let lda=(rmat.len()/ncols) as i32;
                    //let ldb=nrows as i32;
                    let mut info = 0 as i32;

                    /*
                    if check_zero_diagonal(lda as usize,ncols,&rmat){
                        println!("{:?}",*level);
                        println!("{:?}",node);
                        print_mat(lda as usize,ncols,&rmat);
                        panic!("Zero found on diagonal of an R matrix (input sparse matrix was singular)");
                    }
                    */


                    F::xtrtrs(uplo,trans,diag,n,nrhs as i32,&rmat,lda,&mut tmpd,nrowsd as i32,&mut info);



                    if info!=0{
                        print_diag(lda as usize,ncols,&rmat);
                    }
                    assert_eq!(info,0);
                }



                if nrows>0{
                    //Update remaining values
                    //out<----out - A*B
                    {

                        let transa : u8 = b'N';
                        let transb : u8 = b'N';
                        let m = perm.len() as i32;
                        let n = nrhs as i32;
                        let k = self.nodes[node].len() as i32;
                        let amat=&self.triu_num[node];
                        let bmat=&tmpd;

                        //println!("m: {:?},   n:  {:?},    k:   {:?}",m,n,k);
                        //println!("{:?}",*level);
                        //println!("{:?}",node);
                        assert_eq!(bmat.len(),(k*n) as usize);
                        assert_eq!(amat.len(),(m*k) as usize);

                        F::xgemm(transa,transb,m,n,k,F::zero()-F::one(),amat,m,bmat,k,F::one(),&mut tmp,m);
                    }
                }



                //Permute data from temporary array back into in/out array
                permute_in(self.nrows,inout,&perm,&tmp);
                //Permute diagonal block data from temporary array back into in/out array
                permute_in(self.nrows,inout,&permd,&tmpd);
            }
        }

    }
    //total megabytes allocated to this node
    fn node_mb(&self,i : usize)->usize{
        use std::mem::size_of;
        let triu_sz=self.triu_num[i].len()*size_of::<F>();
        let tril_sz=self.tril_num[i].len()*size_of::<F>();
        (triu_sz+tril_sz)/(1024*1024)
    }

    //Plot tree to dotfile using children relationships, annotate with bytes
    //allocated to that node
    pub fn to_graphviz_children_bytes(&self,fname : String) -> std::io::Result<()>{
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(fname)?;
        write!(&mut file,"digraph g {{ \n")?;

        for (i,_c) in self.children.iter().enumerate(){
            let mb = self.node_mb(i);
            write!(&mut file,"  {} [label=\"{}MB\"]\n",i,mb)?;
        }
        for (i,c) in self.children.iter().enumerate(){
            if let Some((c1,c2)) = c{
                write!(&mut file,"  {} -> {}\n",i,c1)?;
                write!(&mut file,"  {} -> {}\n",i,c2)?;
            }
        }
        write!(&mut file,"}} \n")?;
        Ok(())
    }

    //Plot tree to dotfile using children relationships, annotate with time
    //taken to complete all factorization tasks associated with that node
    pub fn to_graphviz_children_times(&self,fname : String) -> std::io::Result<()>{
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(fname)?;
        write!(&mut file,"digraph g {{ \n")?;

        for (i,_c) in self.children.iter().enumerate(){
            let time = self.times[i];
            write!(&mut file,"  {} [label=\"{}s\"]\n",i,time)?;
        }
        for (i,c) in self.children.iter().enumerate(){
            if let Some((c1,c2)) = c{
                write!(&mut file,"  {} -> {}\n",i,c1)?;
                write!(&mut file,"  {} -> {}\n",i,c2)?;
            }
        }
        write!(&mut file,"}} \n")?;
        Ok(())
    }

    //Print aggregate byte consumption in MB by tree level
    pub fn print_level_bytes(&self) -> (){
        for (j,level) in self.levels.iter().enumerate(){
            let mb=level.iter().map(|&i|self.node_mb(i)).fold(0,|acc,x|acc+x);
            println!("level {},   MB:   {}",j,mb);
        }
    }

    //Print aggregate time in seconds by tree level
    pub fn print_level_times(&self) -> (){
        for (j,level) in self.levels.iter().enumerate(){
            let time=level.iter().map(|&i|self.times[i]).fold(0.0,|acc,x|acc+x);
            println!("level {},   MB:   {}",j,time);
        }
    }








    pub fn get_nrows(&self) -> usize{
        self.nrows
    }
    pub fn get_ncols(&self) -> usize{
        self.ncols
    }
    pub fn get_children(&self) -> Vec<Option<(usize,usize)>>{
        self.children.clone()
    }
    pub fn get_parents(&self) -> Vec<Option<usize>>{
        self.parents.clone()
    }

}



#[cfg(test)]
mod tests {
    use crate::sparse::CSCSparse;
    use crate::dtree::DissectionTree;
    use crate::numeric::SparseQR;
    use crate::gallery::laplace2d;
    use crate::gallery::arbstencil2d;
    use crate::numeric::permute_out;
    use crate::numeric::permute_in;


    #[test]
    fn arbstencil2d_solve_f64_verysmall(){
        use num_traits::Zero;
        type F=f64;
        let mx = 16;
        let my = 16;
        let m = mx*my;
        let maxnodes=200;
        //Assemble matrix
        let lap=arbstencil2d::<F>(mx,my,&[-5.0,1.1,-1.2,0.5,3.0]);
        //Get squared graph
        let g = lap.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&lap);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).map(|x| (x as F).sin()+2.0).collect();
        let mut b = vec![F::zero();soln.len()];
        lap.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<F> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).abs()/x.abs()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        assert!(*maxerr < 1e-11);
    }





    #[test]
    fn arbstencil2d_solve_f32(){
        use num_traits::Zero;
        type F=f32;
        let mx = 64;
        let my = 64;
        let m = mx*my;
        let maxnodes=200;
        //Assemble matrix
        let lap=arbstencil2d::<F>(mx,my,&[-5.0,1.1,-1.2,0.5,3.0]);
        //Get squared graph
        let g = lap.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&lap);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).map(|x| (x as F).sin()+2.0).collect();
        let mut b = vec![F::zero();soln.len()];
        lap.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<F> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).abs()/x.abs()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        assert!(*maxerr < 1e-4);
    }

    #[test]
    fn arbstencil2d_solve_f64(){
        use num_traits::Zero;
        type F=f64;
        let mx = 64;
        let my = 64;
        let m = mx*my;
        let maxnodes=200;
        //Assemble matrix
        let lap=arbstencil2d::<F>(mx,my,&[-5.0,1.1,-1.2,0.5,3.0]);
        //Get squared graph
        let g = lap.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&lap);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).map(|x| (x as F).sin()+2.0).collect();
        let mut b = vec![F::zero();soln.len()];
        lap.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<F> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).abs()/x.abs()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        assert!(*maxerr < 1e-11);
    }

    #[test]
    fn arbstencil2d_solve_c32(){
        use num_complex::Complex;
        use num_traits::Zero;
        type R=f32;
        type F=Complex<R>;
        let mx = 64;
        let my = 64;
        let m = mx*my;
        let maxnodes=200;
        //Assemble matrix
        let lap=arbstencil2d::<F>(mx,my,&[F::new(-8.0,3.0),F::new(-1.1,0.3),F::new(0.9,1.2),F::new(1.2,1.1),F::new(0.4,0.3)]);
        //Get squared graph
        let g = lap.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&lap);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).zip(0..m).map(|(xr,xi)| {
            F::new((xr as R).sin(),(xi as R).cos())
        }).collect();

        let mut b = vec![F::zero();soln.len()];
        lap.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<R> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).norm()/x.norm()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        assert!(*maxerr < 1e-4);
    }

    #[test]
    fn arbstencil2d_solve_c64(){
        use num_complex::Complex;
        use num_traits::Zero;
        type R=f64;
        type F=Complex<R>;
        let mx = 64;
        let my = 64;
        let m = mx*my;
        let maxnodes=200;
        //Assemble matrix
        let lap=arbstencil2d::<F>(mx,my,&[F::new(-8.0,3.0),F::new(-1.1,0.3),F::new(0.9,1.2),F::new(1.2,1.1),F::new(0.4,0.3)]);
        //Get squared graph
        let g = lap.to_metis_graph();
        let g2 = g.square();
        //Symbolic factorization
        let dtree = DissectionTree::new(&g2,maxnodes);
        //Numeric factorization, consumes input tree
        let mut fact = SparseQR::<F>::new(dtree,&lap);

        //Make some data to solve with
        let soln : Vec<F> = (0..m).zip(0..m).map(|(xr,xi)| {
            F::new((xr as R).sin(),(xi as R).cos())
        }).collect();

        let mut b = vec![F::zero();soln.len()];
        lap.eval(&soln,&mut b);

        //Now go back and solve
        fact.solve(&mut b);

        //Result should be close to `soln`
        let mut errs : Vec<R> = soln.iter().zip(b.iter()).map(|(&x,&y)|(x-y).norm()/x.norm()).collect();
        errs.sort_by(|x,y|x.partial_cmp(y).unwrap());
        let maxerr = errs.last().unwrap();
        assert!(*maxerr < 1e-12);
    }







    #[test]
    fn permute_out_test(){
        let nrows=32;
        let ncols=32;
        let m=nrows*ncols;
        let vs  : Vec<usize> = (0..m).collect();
        let p : Vec<usize> = vec![1,5,9,10];
        let pvs = permute_out(nrows,&vs,&p);

        assert_eq!(pvs.len(),ncols*p.len());

        for (k,rs) in pvs.chunks_exact(p.len()).into_iter().enumerate(){
            assert_eq!( rs[0],  1 + k*nrows);
            assert_eq!( rs[1],  5 + k*nrows);
            assert_eq!( rs[2],  9 + k*nrows);
            assert_eq!( rs[3], 10 + k*nrows);
        }
    }

    #[test]
    fn permute_in_test(){
        let nrows=10;
        let ncols=32;
        let m=nrows*ncols;
        let vs  : Vec<usize> = (0..m).collect();
        let mut t = vec![0;m];
        let p1 : Vec<usize> = vec![0,2,4,6,8];
        let p2 : Vec<usize> = vec![1,3,5,7,9];

        let t1 = permute_out(nrows,&vs,&p1);
        let t2 = permute_out(nrows,&vs,&p2);

        permute_in(nrows,&mut t,&p1,&t1);
        permute_in(nrows,&mut t,&p2,&t2);

        for (x,y) in t.iter().zip(vs.iter()){
            assert_eq!(*x,*y);
        }
    }





    #[test]
    fn stencil_solve(){
        let mx=16;
        let my=16;
        let m=mx*my;


        let a = laplace2d::<f64>(mx,my);
        //let b = eye::<f64>(m);

        let g = a.to_metis_graph();
        let g2 = g.square();
        let dtree = DissectionTree::new(&g2,32);
        let mut fact = SparseQR::new(dtree,&a);
        //let mut fact = SparseQR::new(dtree,&b);

        let rs : Vec<usize> = (0..m).collect();
        let cs : Vec<usize> = (0..10).collect();
        let mut ys = a.slice_copy(&rs,&cs);

        fact.solve(&mut ys);

        for (k,vals) in ys.chunks_exact(m).into_iter().enumerate(){
            for (j,r) in vals.into_iter().enumerate(){
                if k==j{
                    assert!( (*r-1.0).abs() < 1e-14);
                }
                else{
                    assert!( *r < 1e-14 );
                }
            }
        }
    }



    #[test]
    fn construct_stencil_factorization(){
        let mut offsets = Vec::<usize>::new();
        let mut rids = Vec::<usize>::new();
        let mut vals = Vec::<f64>::new();
        let mut offs=0;

        let mx=32 as usize;
        let my=32 as usize;
        let nrows=mx*my;
        let ncols=mx*my;
        let id = |ix : usize,iy : usize|{iy + my*ix};
        offsets.push(offs);
        for ix in 0..mx{
            for iy in 0..my{
                rids.push(id(ix,iy));
                vals.push(1.0);
                offs+=1;
                if ix>0{
                    rids.push(id(ix-1,iy));
                    vals.push(1.0);
                    offs+=1;
                }
                if ix<mx-1{
                    rids.push(id(ix+1,iy));
                    vals.push(1.0);
                    offs+=1;
                }
                if iy>0{
                    rids.push(id(ix,iy-1));
                    vals.push(1.0);
                    offs+=1;
                }
                if iy<my-1{
                    rids.push(id(ix,iy+1));
                    vals.push(1.0);
                    offs+=1;
                }
                offsets.push(offs);
            }
        }

        let a = CSCSparse::new(nrows,ncols,offsets,rids,vals);
        let g = a.to_metis_graph();
        let g2 = g.square();
        let dtree = DissectionTree::new(&g2,30);
        let _fact = SparseQR::new(dtree,&a);
    }


}
