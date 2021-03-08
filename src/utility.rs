use lapack::{sgeqrf,sormqr,strtrs};
use lapack::{dgeqrf,dormqr,dtrtrs};
use lapack::{cgeqrf,cunmqr,ctrtrs};
use lapack::{zgeqrf,zunmqr,ztrtrs};
use num_complex::Complex;

pub trait Lapack{
    type F;
    //TODO: `a` actually gets written here, but restored after xmqr completes. This reference needs to be a mutable reference.
    fn xmqr(side : u8,trans : u8,m : i32,n : i32,k : i32,a : &[Self::F],lda : i32,tau : &[Self::F],c : &mut [Self::F],ldc : i32,work : &mut [Self::F],lwork : i32,info : &mut i32);
    fn xgeqrf(m: i32,n: i32,a: &mut [Self::F],lda: i32,tau: &mut [Self::F],work: &mut [Self::F],lwork: i32,info: &mut i32);
    fn xtrtrs(uplo: u8,trans: u8,diag: u8,n: i32,nrhs: i32,a: &[Self::F],lda: i32,b: &mut [Self::F],ldb: i32,info: &mut i32);
}



impl Lapack for f32{
    type F=f32;
    fn xmqr(side : u8,trans : u8,m : i32,n : i32,k : i32,a : &[Self::F],lda : i32,tau : &[Self::F],c : &mut [Self::F],ldc : i32,work : &mut [Self::F],lwork : i32,info : &mut i32){
        unsafe{sormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info);}
    }
    fn xgeqrf(m: i32,n: i32,a: &mut [Self::F],lda: i32,tau: &mut [Self::F],work: &mut [Self::F],lwork: i32,info: &mut i32){
        unsafe{sgeqrf(m,n,a,lda,tau,work,lwork,info);}
    }
    fn xtrtrs(uplo: u8,trans: u8,diag: u8,n: i32,nrhs: i32,a: &[Self::F],lda: i32,b: &mut [Self::F],ldb: i32,info: &mut i32){
        unsafe{strtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info);}
    }
}

impl Lapack for f64{
    type F=f64;
    fn xmqr(side : u8,trans : u8,m : i32,n : i32,k : i32,a : &[Self::F],lda : i32,tau : &[Self::F],c : &mut [Self::F],ldc : i32,work : &mut [Self::F],lwork : i32,info : &mut i32){
        unsafe{dormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info);}
    }
    fn xgeqrf(m: i32,n: i32,a: &mut [Self::F],lda: i32,tau: &mut [Self::F],work: &mut [Self::F],lwork: i32,info: &mut i32){
        unsafe{dgeqrf(m,n,a,lda,tau,work,lwork,info);}
    }
    fn xtrtrs(uplo: u8,trans: u8,diag: u8,n: i32,nrhs: i32,a: &[Self::F],lda: i32,b: &mut [Self::F],ldb: i32,info: &mut i32){
        unsafe{dtrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info);}
    }
}

impl Lapack for Complex<f32>{
    type F=Complex<f32>;
    fn xmqr(side : u8,trans : u8,m : i32,n : i32,k : i32,a : &[Self::F],lda : i32,tau : &[Self::F],c : &mut [Self::F],ldc : i32,work : &mut [Self::F],lwork : i32,info : &mut i32){
        let ctrans = if trans==b'T' {b'C'} else {trans};
        unsafe{cunmqr(side,ctrans,m,n,k,a,lda,tau,c,ldc,work,lwork,info);}
    }
    fn xgeqrf(m: i32,n: i32,a: &mut [Self::F],lda: i32,tau: &mut [Self::F],work: &mut [Self::F],lwork: i32,info: &mut i32){
        unsafe{cgeqrf(m,n,a,lda,tau,work,lwork,info);}
    }
    fn xtrtrs(uplo: u8,trans: u8,diag: u8,n: i32,nrhs: i32,a: &[Self::F],lda: i32,b: &mut [Self::F],ldb: i32,info: &mut i32){
        let ctrans = if trans==b'T' {b'C'} else {trans};
        unsafe{ctrtrs(uplo,ctrans,diag,n,nrhs,a,lda,b,ldb,info);}
    }
}

impl Lapack for Complex<f64>{
    type F=Complex<f64>;
    fn xmqr(side : u8,trans : u8,m : i32,n : i32,k : i32,a : &[Self::F],lda : i32,tau : &[Self::F],c : &mut [Self::F],ldc : i32,work : &mut [Self::F],lwork : i32,info : &mut i32){
        let ctrans = if trans==b'T' {b'C'} else {trans};
        unsafe{zunmqr(side,ctrans,m,n,k,a,lda,tau,c,ldc,work,lwork,info);}
    }
    fn xgeqrf(m: i32,n: i32,a: &mut [Self::F],lda: i32,tau: &mut [Self::F],work: &mut [Self::F],lwork: i32,info: &mut i32){
        unsafe{zgeqrf(m,n,a,lda,tau,work,lwork,info);}
    }
    fn xtrtrs(uplo: u8,trans: u8,diag: u8,n: i32,nrhs: i32,a: &[Self::F],lda: i32,b: &mut [Self::F],ldb: i32,info: &mut i32){
        let ctrans = if trans==b'T' {b'C'} else {trans};
        unsafe{ztrtrs(uplo,ctrans,diag,n,nrhs,a,lda,b,ldb,info);}
    }
}
