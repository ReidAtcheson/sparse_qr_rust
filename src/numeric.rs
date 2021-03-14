use crate::dtree::DissectionTree;
use crate::sparse::CSCSparse;
use crate::utility::Lapack;
use std::collections::BTreeMap;
use num_traits::Num;
use num_traits::cast::ToPrimitive;


fn check_zero_diagonal<F : Num>(nrows : usize, ncols : usize, a : &Vec<F>) -> bool{
    let mut out : bool = false;
    for i in 0..ncols{
        let j=i;
        if a[i+nrows*j] == F::zero() {
            out=true;
        }
    }
    out
}


fn panic_zero_diagonal<F : Num>(nrows : usize, ncols : usize, a : &Vec<F>) -> (){
    for i in 0..ncols{
        let j=i;
        assert!( a[i+nrows*j] != F::zero() );
    }
}

fn print_mat<F : std::fmt::Display>(nrows : usize, ncols : usize,a : &Vec<F>)->(){
    for i in 0..nrows{
        for j in 0..ncols{
            print!("{} ",a[i+nrows*j]);
        }
        println!();
    }
}
fn print_diag<F : std::fmt::Display>(nrows : usize, ncols : usize,a : &Vec<F>)->(){
    for i in 0..nrows{
        let j=i;
        print!("{}\n",a[i+nrows*j]);
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
} 


impl <F : Lapack<F=F>+Num+ToPrimitive+Copy + std::fmt::Display> SparseQR<F>{
    pub fn new(dtree : DissectionTree, mat : &CSCSparse<F>) -> Self{

        let nnodes=dtree.nodes.len();
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
        for mut tri in tril.iter_mut(){
            tri.sort();
        }
        for mut tri in triu.iter_mut(){
            tri.sort();
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
                //Gather lower triangular blocks into a stacked matrix
                let ncols=dtree.nodes[node].len();
                let nrows=tril[node].iter().map(|&k|dtree.nodes[k].len()).fold(0,|acc,x|acc+x);
                tril_num[node]=vec![F::zero();nrows*ncols];
                let trilmat=&mut tril_num[node];
                let nmap=&temp[node];
                let mut offs=0;
                for tri in tril[node].iter().cloned(){
                    let lnrows=dtree.nodes[tri].len();
                    let m = nmap.get(&tri).unwrap();
                    /*
                    println!("==================");
                    print_mat(lnrows,ncols,m);
                    println!("==================");
                    */
                    for (xs,ys) in trilmat.chunks_exact_mut(nrows).zip(m.chunks_exact(lnrows)){
                        for (x,y) in (&mut xs[offs..offs+lnrows]).into_iter().zip(ys){
                            *x=*y;
                        }
                    }
                    offs+=lnrows;
                    /*
                    print_mat(nrows,ncols,trilmat);
                    */
                }


                //First call queries optimal workspace size
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
                        if check_zero_diagonal(nrows,ncols,&trilmat){
                            println!("{:?}",*level);
                            println!("{:?}",node);
                            print_mat(nrows,ncols,&trilmat);
                            panic!("Zero found on diagonal of an R matrix (input sparse matrix was singular)");
                        }
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
                    let mut triumat = vec![F::zero();triu_ncols*triu_nrows];
                    {
                        let mut offs=0;
                        for k in tril[node].iter(){
                            let tmat=temp[p].get(k).unwrap();
                            let lnrows=dtree.nodes[*k].len();
                            for (xs,ys) in triumat.chunks_exact_mut(nrows).zip(tmat.chunks_exact(lnrows)){
                                for (x,y) in (&mut xs[offs..offs+lnrows]).into_iter().zip(ys){
                                    *x=*y;
                                }
                            }
                            offs+=lnrows;
                        }
                    }


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
                        F::xmqr(side,trans,triu_nrows as i32,triu_ncols as i32,tau.len() as i32,trilmat.as_slice(),nrows as i32,tau.as_slice(),
                        triumat.as_mut_slice(),triu_nrows as i32,work.as_mut_slice(),lwork,&mut info);
                        assert_eq!(info,0);
                    };

                    //Unpack matrices
                    {

                        let mut offs=0;
                        for k in tril[node].iter(){
                            let tmat=temp[p].get_mut(k).unwrap();
                            let lnrows=dtree.nodes[*k].len();
                            for (xs,ys) in triumat.chunks_exact(nrows).zip(tmat.chunks_exact_mut(lnrows)){
                                for (x,y) in xs[offs..offs+lnrows].into_iter().zip(ys.iter_mut()){
                                    *y=*x;
                                }
                            }
                            offs+=lnrows;
                        }
                    }
                    mp=dtree.parents[p];
                }

                //Save tau scalars for the houesholder reflectors
                taus[node]=tau;
            }
        }

        //Factorization complete but now we need to re-pack the upper
        //triangular matrices
        let mut triu_num  = vec![Vec::<F>::new();nnodes];
        for (k,tri) in triu.iter().enumerate(){
            let ncols=dtree.nodes[k].len();
            let nrows=tri.iter().map(|&k|dtree.nodes[k].len()).fold(0,|acc,x|acc+x);
            let mut packed = vec![F::zero();nrows*ncols];
            let mut offs=0;
            for t in tri.iter(){
                let lnrows=dtree.nodes[*t].len();
                let tmat=temp[k].get(t).unwrap();
                for (xs,ys) in packed.chunks_exact_mut(nrows).zip(tmat.chunks_exact(lnrows)){
                    for (x,y) in (&mut xs[offs..offs+lnrows]).into_iter().zip(ys.iter()){
                        *x=*y;
                    }
                }
                offs+=lnrows;
            }
        }


        SparseQR { parents : dtree.parents, children : dtree.children, levels : dtree.levels, nodes : dtree.nodes,
        nrows : mat.get_nrows(), ncols : mat.get_ncols(), triu : triu, tril : tril, triu_num : triu_num, tril_num : tril_num,taus : taus
        }
    }

    pub fn solve(&mut self,inout : &mut [F]) -> (){
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
                let ncols=self.nodes[node].len();
                //Permute active data into temporary array
                let mut tmp = {
                    let mut tmp = vec![F::zero();nrows*nrhs];
                    //The temporary array and actual input should have same number of columns
                    assert_eq!(tmp.len()/nrows,inout.len()/self.nrows);
                    //Loop over columns
                    for (xs,ys) in tmp.as_mut_slice().chunks_exact_mut(nrows).zip(inout.chunks_exact(self.nrows)){
                        //Loop over entries of column
                        for (i,x) in xs.iter_mut().enumerate(){
                            let pi=perm[i];
                            *x=ys[pi];
                        }
                    }
                    tmp
                };
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
                for (xs,ys) in tmp.as_slice().chunks_exact(nrows).zip(inout.chunks_exact_mut(self.nrows)){
                    for (i,x) in xs.iter().enumerate(){
                        let pi=perm[i];
                        ys[pi]=*x;
                    }
                }
            }
        }
        //Now `inout` has been overwritten with ` Q^T * inout`. We now do back-substitution
        //which is similar to previous loop but in reversed order starting
        //with the top-most node of the nested dissection tree
        for level in self.levels.iter(){
            for node in level.iter().cloned(){
                //Get permutation vector corresponding to upper triangular part
                let perm : Vec<usize> = self.tril[node].iter().map(|&n|self.nodes[n].clone()).flatten().collect();
                //Permutation corresponding to diagonal block
                let permd : Vec<usize> = self.nodes[node].clone();
                let nrows=perm.len();
                let nrowsd=permd.len();
                let ncols=self.nodes[node].len();
                //Permute upper triangular active data into temporary array
                let mut tmp = {
                    let mut tmp = vec![F::zero();nrows*nrhs];
                    //The temporary array and actual input should have same number of columns
                    assert_eq!(tmp.len()/nrows,inout.len()/self.nrows);
                    //Loop over columns
                    for (xs,ys) in tmp.as_mut_slice().chunks_exact_mut(nrowsd).zip(inout.chunks_exact(self.nrows)){
                        //Loop over entries of column
                        for (i,x) in xs.iter_mut().enumerate(){
                            let pi=perm[i];
                            *x=ys[pi];
                        }
                    }
                    tmp
                };

                //Permute diagonal block active data into temporary array
                let mut tmpd = {
                    let mut tmpd = vec![F::zero();nrowsd*nrhs];
                    //The temporary array and actual input should have same number of columns
                    assert_eq!(tmpd.len()/nrowsd,inout.len()/self.nrows);
                    //Loop over columns
                    for (xs,ys) in tmpd.as_mut_slice().chunks_exact_mut(nrowsd).zip(inout.chunks_exact(self.nrows)){
                        //Loop over entries of column
                        for (i,x) in xs.iter_mut().enumerate(){
                            let pi=permd[i];
                            *x=ys[pi];
                        }
                    }
                    tmpd
                };


                //Back-substitution of diagonal block
                {
                    let rmat = &self.tril_num[node];
                    let uplo : u8 =b'U';
                    let trans : u8 =b'N';
                    let diag : u8 =b'N';
                    let n=ncols as i32;
                    assert_eq!(rmat.len() % ncols, 0);
                    let lda=(rmat.len()/ncols) as i32;
                    let ldb=nrows as i32;
                    let mut info = 0 as i32;

                    if check_zero_diagonal(lda as usize,ncols,&rmat){
                        println!("{:?}",*level);
                        println!("{:?}",node);
                        print_mat(lda as usize,ncols,&rmat);
                        panic!("Zero found on diagonal of an R matrix (input sparse matrix was singular)");
                    }


                    F::xtrtrs(uplo,trans,diag,n,nrhs as i32,&rmat,lda,&mut tmpd,nrowsd as i32,&mut info);



                    if info!=0{
                        print_diag(lda as usize,ncols,&rmat);
                    }
                    assert_eq!(info,0);
                }


                //Update remaining values



                //Permute data from temporary array back into in/out array
                for (xs,ys) in tmp.as_slice().chunks_exact(nrows).zip(inout.chunks_exact_mut(self.nrows)){
                    for (i,x) in xs.iter().enumerate(){
                        let pi=perm[i];
                        ys[pi]=*x;
                    }
                }
                //Permute diagonal block data from temporary array back into in/out array
                for (xs,ys) in tmpd.as_slice().chunks_exact(nrowsd).zip(inout.chunks_exact_mut(self.nrows)){
                    //Loop over entries of column
                    for (i,x) in xs.iter().enumerate(){
                        let pi=permd[i];
                        ys[pi]=*x;
                    }
                }





            }
        }


    }
}



#[cfg(test)]
mod tests {
    use crate::sparse::CSCSparse;
    use crate::metis::MetisGraph;
    use crate::dtree::DissectionTree;
    use crate::numeric::SparseQR;
    use crate::gallery::laplace2d;


    #[test]
    fn stencil_solve(){
        let mx=16;
        let my=16;
        let m=mx*my;


        let a = laplace2d::<f64>(mx,my);
        let g = a.to_metis_graph();
        let g2 = g.square();
        let dtree = DissectionTree::new(&g2,32);
        let mut fact = SparseQR::new(dtree,&a);

        let rs : Vec<usize> = (0..m).collect();
        let cs : Vec<usize> = (0..10).collect();
        let mut ys = a.slice_copy(&rs,&cs);
        fact.solve(&mut ys);
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
        let fact = SparseQR::new(dtree,&a);
    }


}
