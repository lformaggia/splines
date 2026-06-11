/**
 * @file concepts.hpp
 * @brief C++20 concepts that constrain numeric template parameters.
 *
 * This header defines the `RealType` concept used throughout the library
 * to restrict template type parameters to IEEE 754 floating-point types
 * (`float`, `double`, `long double`).
 */
#pragma once
#include <concepts>
#include <ranges>

namespace splines {

/**
 * @brief Constrains a type to a floating-point (real) scalar.
 *
 * Satisfied by `float`, `double`, and `long double`.  Every library template
 * that accepts a scalar type uses this concept to produce a clear compilation
 * error when an unsupported type is supplied.
 *
 * @tparam T  The type to test.
 *
 * @par Example
 * @code
 * template<RealType Real>
 * Real squareRoot(Real x) { return std::sqrt(x); }
 * @endcode
 */
template<typename T>
concept RealType = std::floating_point<T>;

} // namespace splines
