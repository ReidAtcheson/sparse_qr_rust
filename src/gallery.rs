use num_traits::Num;
use num_traits::cast::NumCast;

use crate::sparse::CSCSparse;


pub fn eye<F : Num+NumCast+Copy>(m : usize) -> CSCSparse<F> {
    let nrows=m;
    let ncols=m;
    let mut offsets = Vec::<usize>::new();
    let mut rids = Vec::<usize>::new();
    let mut vals = Vec::<F>::new();
    let mut offs=0;
    offsets.push(offs);
    for i in 0..m{
        rids.push(i);
        vals.push(F::from(1.0).unwrap());
        offs+=1;
        offsets.push(offs);
    }
    CSCSparse::<F>::new(nrows,ncols,offsets,rids,vals)
}

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

pub fn arbstencil2d<F : Num+NumCast+Copy>(mx : usize,my : usize,cs : &[F;5]) -> CSCSparse<F> {
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
            vals.push(cs[0]);
            offs+=1;
            if ix>0{
                rids.push(id(ix-1,iy));
                vals.push(cs[1]);
                offs+=1;
            }
            if ix<mx-1{
                rids.push(id(ix+1,iy));
                vals.push(cs[2]);
                offs+=1;
            }
            if iy>0{
                rids.push(id(ix,iy-1));
                vals.push(cs[3]);
                offs+=1;
            }
            if iy<my-1{
                rids.push(id(ix,iy+1));
                vals.push(cs[4]);
                offs+=1;
            }
            offsets.push(offs);
        }
    }
    CSCSparse::<F>::new(nrows,ncols,offsets,rids,vals)
}

pub fn arbstencil3d<F : Num+NumCast+Copy>(mx : usize,my : usize,mz : usize,cs : &[F;7]) -> CSCSparse<F> {
    let nrows=mx*my*mz;
    let ncols=mx*my*mz;
    let mut offsets = Vec::<usize>::new();
    let mut rids = Vec::<usize>::new();
    let mut vals = Vec::<F>::new();
    let mut offs=0;
    offsets.push(offs);
    let id = |ix : usize,iy : usize,iz : usize|{iz + mz*(iy + my*ix)};
    for ix in 0..mx{
        for iy in 0..my{
            for iz in 0..mz{
                rids.push(id(ix,iy,iz));
                vals.push(cs[0]);
                offs+=1;
                if ix>0{
                    rids.push(id(ix-1,iy,iz));
                    vals.push(cs[1]);
                    offs+=1;
                }
                if ix<mx-1{
                    rids.push(id(ix+1,iy,iz));
                    vals.push(cs[2]);
                    offs+=1;
                }
                if iy>0{
                    rids.push(id(ix,iy-1,iz));
                    vals.push(cs[3]);
                    offs+=1;
                }
                if iy<my-1{
                    rids.push(id(ix,iy+1,iz));
                    vals.push(cs[4]);
                    offs+=1;
                }
                if iz>0{
                    rids.push(id(ix,iy,iz-1));
                    vals.push(cs[5]);
                    offs+=1;
                }
                if iz<mz-1{
                    rids.push(id(ix,iy,iz+1));
                    vals.push(cs[6]);
                    offs+=1;
                }



                offsets.push(offs);
            }
        }
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

pub fn laplace3d<F : Num+NumCast+Copy>(mx : usize,my : usize,mz : usize) -> CSCSparse<F> {
    let nrows=mx*my*mz;
    let ncols=mx*my*mz;
    let mut offsets = Vec::<usize>::new();
    let mut rids = Vec::<usize>::new();
    let mut vals = Vec::<F>::new();
    let mut offs=0;
    offsets.push(offs);
    let id = |ix : usize,iy : usize,iz : usize|{iz + mz*(iy + my*ix)};
    for ix in 0..mx{
        for iy in 0..my{
            for iz in 0..mz{
                rids.push(id(ix,iy,iz));
                vals.push(F::from(6.0).unwrap());
                offs+=1;
                if ix>0{
                    rids.push(id(ix-1,iy,iz));
                    vals.push(F::from(-1.0).unwrap());
                    offs+=1;
                }
                if ix<mx-1{
                    rids.push(id(ix+1,iy,iz));
                    vals.push(F::from(-1.0).unwrap());
                    offs+=1;
                }
                if iy>0{
                    rids.push(id(ix,iy-1,iz));
                    vals.push(F::from(-1.0).unwrap());
                    offs+=1;
                }
                if iy<my-1{
                    rids.push(id(ix,iy+1,iz));
                    vals.push(F::from(-1.0).unwrap());
                    offs+=1;
                }
                if iz>0{
                    rids.push(id(ix,iy,iz-1));
                    vals.push(F::from(-1.0).unwrap());
                    offs+=1;
                }
                if iz<mz-1{
                    rids.push(id(ix,iy,iz+1));
                    vals.push(F::from(-1.0).unwrap());
                    offs+=1;
                }
                offsets.push(offs);
            }
        }
    }
    CSCSparse::<F>::new(nrows,ncols,offsets,rids,vals)
}



#[cfg(test)]
mod tests {
    use crate::gallery::{laplace2d,eye};


    #[test]
    fn gallery_id(){
        let mx=16;
        let my=16;
        let m=mx*my;
        let a = eye::<f64>(m);
        let b = laplace2d::<f64>(mx,my);

        let bvec = b.slice_copy(&(0..m).collect(),&(0..10).collect());
        let mut out = vec![0.0;m*10];

        a.eval(&bvec,&mut out);

        for (x,y) in bvec.iter().zip(out.iter()){
            assert_eq!(*x,*y);
        }
    }

}

