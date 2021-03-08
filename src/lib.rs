extern crate blas;
extern crate openblas_src;

pub mod metis;
pub mod dtree;
pub mod sparse;
pub mod numeric;
pub mod utility;
pub mod gallery;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
