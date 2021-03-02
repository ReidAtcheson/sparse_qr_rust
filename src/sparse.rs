use num_traits::Num;

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
        let mut offs = offsets;
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
    pub fn eval(&self,x : &[F], y : &mut [F]) -> (){
        assert!(x.len()>0);
        assert!(y.len()>0);
        assert!(x.len() % self.nrows == 0);
        assert!(y.len() % self.nrows == 0);
        assert!(x.len() == y.len());
        for yv in y.iter_mut(){
            *yv=F::zero();
        }

        for xc in x.chunks_exact(self.nrows){
            for j in 0..xc.len(){
                let beg=self.offsets[j];
                let end=self.offsets[j+1];
                for (r,nz) in self.rids[beg..end].iter().zip(self.vals[beg..end].iter()){
                    y[*r] = y[*r]+x[j]*(*nz);
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

}


