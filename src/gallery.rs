use num_traits::Num;
use num_traits::cast::NumCast;

use crate::sparse::CSCSparse;


pub fn laplace1d<F : Num+NumCast+Copy>(m : usize) -> CSCSparse<F> {
    let nrows=m;
    let ncols=m;
    let mut offsets = Vec::<usize>::new();
    let mut rids = Vec::<usize>::new();
    let mut vals = Vec::<F>::new();

    let mut offs=0;
    offsets.push(offs);
    for i in 0..m{
        if i>0{
            rids.push(i-1);
            vals.push(F::from(-1.0).unwrap());
            offs+=1;
        }
        if i<m-1{
            rids.push(i+1);
            vals.push(F::from(-1.0).unwrap());
            offs+=1;
        }
        rids.push(i);
        vals.push(F::from(2.0).unwrap());
        offs+=1;
        offsets.push(offs);
    }
    CSCSparse::<F>::new(nrows,ncols,offsets,rids,vals)
}

pub fn laplace2d<F : Num+NumCast+Copy>(mx : usize,my : usize) -> CSCSparse<F> {
    let nrows=mx*my;
    let ncols=mx*my;
    let mut offsets = Vec::<usize>::new();
    let mut rids = Vec::<usize>::new();
    let mut vals = Vec::<F>::new();
    let mut offs=0;
    offsets.push(offs);
    let id = |ix : usize,iy : usize|{iy + my*ix};
    for ix in 0..mx{
        for iy in 0..my{
            rids.push(id(ix,iy));
            vals.push(F::from(4.0).unwrap());
            offs+=1;
            if ix>0{
                rids.push(id(ix-1,iy));
                vals.push(F::from(-1.0).unwrap());
                offs+=1;
            }
            if ix<mx-1{
                rids.push(id(ix+1,iy));
                vals.push(F::from(-1.0).unwrap());
                offs+=1;
            }
            if iy>0{
                rids.push(id(ix,iy-1));
                vals.push(F::from(-1.0).unwrap());
                offs+=1;
            }
            if iy<my-1{
                rids.push(id(ix,iy+1));
                vals.push(F::from(-1.0).unwrap());
                offs+=1;
            }
            offsets.push(offs);
        }
    }
    CSCSparse::<F>::new(nrows,ncols,offsets,rids,vals)
}
