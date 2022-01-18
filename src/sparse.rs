use num_traits::Num;
use crate::metis::MetisGraph;

///A simple compressed column sparse matrix
pub struct CSCSparse<F>{
    nrows : usize,
    ncols : usize,
    //Offsets into column data
    offsets : Vec<usize>,
    //Nonzero row entries
    rids : Vec<usize>,
    //Nonzero values
    vals : Vec<F>
}


impl <F : Num+Copy> CSCSparse<F>{ 

    pub fn panic_if_invalid(&self){
        assert!(self.nrows>0);
        assert!(self.ncols>0);
        //TODO: Replace with "is_sorted" when it merges
        for i in 1..self.offsets.len(){
            assert!(self.offsets[i-1]<=self.offsets[i]);
        }
        for r in self.rids.iter(){
            assert!(*r < self.nrows);
        }
        assert_eq!(self.offsets.len()-1,self.ncols);
        for i in 1..self.offsets.len(){
            let beg=self.offsets[i-1];
            let end=self.offsets[i];
            //Require rows to be sorted in the matrix
            for j in beg+1..end{
                assert!(self.rids[j-1]<self.rids[j]);
            }
        }
        assert_eq!(self.rids.len(),self.vals.len());
    }

    pub fn new(nrows : usize, ncols : usize, offsets : Vec<usize>,rids : Vec<usize>,vals : Vec<F>) -> Self{

        //Put the row indices in sorted order
        let offs = offsets;
        let mut rs = rids;
        let mut vs = vals;
        for i in 1..offs.len(){
            let beg=offs[i-1];
            let end=offs[i];
            let mut perm : Vec<usize> = (0..end-beg).collect();
            perm.sort_by_key(|j|{
                let locrs=&rs[beg..end];
                locrs[*j]
            });
            let mut new_rs : Vec<usize> = vec![0;end-beg];
            let mut new_vs : Vec<F> = vec![F::zero();end-beg];

            for j in 0..end-beg{
                let locrs = &rs[beg..end];
                let locvs = &vs[beg..end];
                new_rs[j] = locrs[perm[j]];
                new_vs[j] = locvs[perm[j]];
            }

            for j in 0..end-beg{
                let locrs = &mut rs[beg..end];
                let locvs = &mut vs[beg..end];
                locrs[j]=new_rs[j];
                locvs[j]=new_vs[j];
            }
        }

        let out = CSCSparse { nrows : nrows, ncols : ncols, offsets : offs,rids : rs, vals : vs};
        out.panic_if_invalid();
        out
    }

    pub fn get_nrows(&self) -> usize{ self.nrows }
    pub fn get_ncols(&self) -> usize{ self.ncols }
    pub fn to_metis_graph(&self) -> MetisGraph{
        let mut offs : i64 = 0;
        let mut xadj = Vec::<i64>::new();
        xadj.push(0);
        let mut adjncy = Vec::<i64>::new();

        for j in 1..self.offsets.len(){
            let beg=self.offsets[j-1];
            let end=self.offsets[j];
            for e in &self.rids[beg..end]{
                if *e != (j-1){
                    adjncy.push(*e as i64);
                    offs+=1;
                }
            }
            xadj.push(offs);
        }
        MetisGraph::new(xadj,adjncy)
    }
    pub fn eval(&self,x : &[F], y : &mut [F]) -> (){
        assert!(x.len()>0);
        assert!(y.len()>0);
        assert!(x.len() % self.nrows == 0);
        assert!(y.len() % self.nrows == 0);
        assert!(x.len() == y.len());
        for yv in y.iter_mut(){
            *yv=F::zero();
        }

        for (xc,yc) in x.chunks_exact(self.nrows).zip(y.chunks_exact_mut(self.nrows)){
            for j in 0..xc.len(){
                let beg=self.offsets[j];
                let end=self.offsets[j+1];
                for (r,nz) in self.rids[beg..end].iter().zip(self.vals[beg..end].iter()){
                    yc[*r]=yc[*r]+xc[j]*(*nz);
                }
            }
        }

    }

    pub fn slice_copy(&self,rows : &Vec<usize>,cols : &Vec<usize>) -> Vec<F>{
        let nrows=rows.len();
        let ncols=cols.len();
        assert!(nrows>0);
        assert!(ncols>0);
        assert!(nrows<=self.nrows);
        assert!(ncols<=self.ncols);
        for r in rows.iter(){
            assert!(*r<self.nrows);
        }
        for c in cols.iter(){
            assert!(*c<self.ncols);
        }

        let mut out = vec![F::zero();nrows*ncols];
        for (j,c) in cols.iter().enumerate(){
            let beg=self.offsets[*c];
            let end=self.offsets[*c+1];
            let rids=&self.rids[beg..end];
            let vals=&self.vals[beg..end];
            for (i,r) in rows.iter().enumerate(){
                if let Ok(k) = rids.binary_search(r){
                    let id=i+nrows*j;
                    out[id]=vals[k];
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use crate::sparse::CSCSparse;
    use crate::gallery::{laplace1d,laplace2d,eye};


    #[test]
    fn sparse_test_slice_id(){
        let mx=16;
        let my=16;
        let m=mx*my;
        let a = eye::<f64>(m);
        let perm=vec![10, 11, 25, 26, 40, 41, 55, 56, 70, 71, 87, 88, 102, 103, 118, 119, 134, 135, 151, 152, 167, 168, 183, 184, 199, 200, 214, 215, 231, 232, 246, 247];
        let asep = a.slice_copy(&perm,&perm);
        for i in 0..perm.len(){
            let j=i;
            assert_eq!(asep[i+perm.len()*j],1.0);
        }
    }

    #[test]
    fn sparse_test_slice1(){
        let m=32;
        let a=laplace1d::<f64>(m);
        let rs : Vec<usize> = (0..5).collect();
        let cs : Vec<usize> = (10..m).collect();
        let b=a.slice_copy(&rs,&cs);
        assert_eq!(b.len(),rs.len()*cs.len());

        for x in b.iter(){
            assert_eq!(*x,0.0);
        }
    }

    #[test]
    fn sparse_test_slice2(){
        let m=32;
        let a=laplace1d::<f64>(m);
        let rs : Vec<usize> = (1..4).collect();
        let cs : Vec<usize> = (2..4).collect();
        let b=a.slice_copy(&rs,&cs);
        assert_eq!(b.len(),rs.len()*cs.len());
        assert_eq!(b[0],-1.0);
    }
    #[test]
    fn sparse_test_slice3(){
        let mx=8;
        let my=8;

        let a=laplace2d::<f64>(mx,my);
        let id = |ix,iy|{iy+my*ix};
        let rs : Vec<usize> = vec![id(3,3)];
        let cs : Vec<usize> = vec![id(3,2),id(3,3),id(3,4),id(2,3),id(4,3)];
        let b=a.slice_copy(&rs,&cs);
        assert_eq!(b.len(),rs.len()*cs.len());

        assert_eq!(b[0],-1.0);
        assert_eq!(b[1],4.0);
        assert_eq!(b[2],-1.0);
        assert_eq!(b[3],-1.0);
        assert_eq!(b[4],-1.0);
    }





    #[test]
    fn sparse_check_rows_sorted(){
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

        let _a = CSCSparse::new(nrows,ncols,offsets,rids,vals);
    }

    #[test]
    fn sparse_make_metis_graph(){
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
        let _g = a.to_metis_graph();
    }


    #[test]
    fn sparse_eval_slice_same(){
        let mut offsets = Vec::<usize>::new();
        let mut rids = Vec::<usize>::new();
        let mut vals = Vec::<f64>::new();
        let mut offs=0;

        let mx=4 as usize;
        let my=4 as usize;
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
                    vals.push(2.0);
                    offs+=1;
                }
                if ix<mx-1{
                    rids.push(id(ix+1,iy));
                    vals.push(3.0);
                    offs+=1;
                }
                if iy>0{
                    rids.push(id(ix,iy-1));
                    vals.push(4.0);
                    offs+=1;
                }
                if iy<my-1{
                    rids.push(id(ix,iy+1));
                    vals.push(5.0);
                    offs+=1;
                }
                offsets.push(offs);
            }
        }

        let a = CSCSparse::new(nrows,ncols,offsets,rids,vals);

        let rows : Vec<usize> = (0..nrows).collect();
        let cols : Vec<usize> = (0..ncols).collect();

        let m1 = a.slice_copy(&rows,&cols);

        let mut m2 : Vec<f64> = vec![0.0;nrows*ncols];
        let mut m3 : Vec<f64> = vec![0.0;nrows*ncols];
        for (i,v) in m2.chunks_exact_mut(nrows).enumerate(){
            v[i]=1.0;
        }
        a.eval(&m2,&mut m3);




        for (x,y) in m1.iter().zip(m3.iter()){
            assert_eq!(*x,*y);
        }

    }


}


