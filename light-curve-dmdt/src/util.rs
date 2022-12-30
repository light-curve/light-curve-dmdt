use itertools::Itertools;

/// Checks if slice is sorted and have no duplicates
///
/// It could be replace with the corresponding standard library method when it will be stabilized
/// <https://github.com/rust-lang/rust/issues/53485>
pub fn is_sorted<T>(a: &[T]) -> bool
where
    T: PartialOrd,
{
    a.iter().tuple_windows().all(|(a, b)| a < b)
}
