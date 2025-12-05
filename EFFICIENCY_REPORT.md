# Efficiency Improvements Report for magenta

This report documents several areas in the codebase where efficiency could be improved.

## Profiling Results

Profiling a simulation with 10,000 people for 10 years with 30-day updates revealed that **98.5% of execution time is spent in C++ code** (via `.Call` / Rcpp), specifically in `Simulation_Update_cpp`. The breakdown shows:

| Component | Time (%) |
|-----------|----------|
| `.Call` (C++ code) | 98.48% |
| `Simulation_Update_cpp` | 97.52% |
| `gc` (garbage collection) | 0.34% |
| R-side functions (`ztrgeomintp`, `ztrnbinom`, etc.) | < 0.5% |

This means that **R-level optimizations have minimal impact** on overall performance. The real bottlenecks are in the C++ simulation loop.

## C++ Optimizations Implemented

### 1. Batch Insert in Multinomial Sampling (probability.cpp)

**Location:** `src/probability.cpp`, `rmultinomN()` function

**Issue:** The original implementation used individual `emplace_back()` calls in a while loop for each draw:
```cpp
while (draw > 0) {
  output.emplace_back(k);
  draw--;
}
```

**Fix:** Replaced with batch `insert()` operations:
```cpp
if (draw > 0) {
  output.insert(output.end(), draw, k);
}
```

**Impact:** Reduces function call overhead and improves cache efficiency when `draw` is large.

### 2. Simplified Biting Frequency Calculation (main_update.cpp)

**Location:** `src/main_update.cpp`, daily biting frequency calculation

**Issue:** The original code used three passes over a vector (generate + transform + adjacent_difference) to compute what is essentially a constant step value:
```cpp
std::generate(temp_biting_frequency_vector.begin(), ...);
std::transform(temp_biting_frequency_vector.begin(), ...);
std::adjacent_difference(temp_biting_frequency_vector.begin(), ...);
```

**Fix:** Replaced with a single `std::fill()` operation:
```cpp
const double biting_step = 1.0 / mosquito_biting_rates[intervention_counter];
std::fill(u_ptr->parameters.g_mosquito_next_biting_day_vector.begin(),
          u_ptr->parameters.g_mosquito_next_biting_day_vector.end(),
          static_cast<int>(biting_step));
```

**Impact:** Reduces three vector passes to one, eliminating unnecessary intermediate computations.

## R-Level Issues (Lower Priority)

The following R-level issues have minimal impact on overall performance but are documented for completeness:

### 1. Inefficient Vector Growth in Random Number Generation Functions

**Location:** `R/stat_genetics_utils.R`, lines 86-116

**Functions:** `ztrgeomintp()` and `ztrnbinom()`

**Issue:** Both functions use while loops with repeated vector concatenation using `c()`. In R, this pattern is inefficient because each concatenation creates a new vector and copies all existing elements.

```r
# Current inefficient pattern:
while(length(ngs) < n) {
  ng2 <- rgeom(n*1.5,prob = 1/mean)
  ngs <- c(ngs,round(ng2*p))  # Creates new vector each iteration
  ngs <- ngs[ngs>0]
}
```

**Impact:** These functions are called in the main simulation pipeline (`pipeline.R` lines 284-285, 390, 599-600) during spatial list creation, which happens at initialization and during annual updates.

**Recommended Fix:** Pre-allocate a larger vector upfront and fill it, or use a more efficient approach that generates enough samples in fewer iterations.

## 2. Deprecated `class()` Usage for Type Checking

**Location:** Multiple files throughout the codebase

**Files affected:**
- `R/mainFunctions.R` (lines 32, 101, 153, 180, 211)
- `R/pipeline_utils.R` (lines 400, 415)
- `R/cluster_utils.R` (lines 14)

**Issue:** The code uses `class(x) != "type"` pattern for type checking, which is deprecated. R documentation recommends using `inherits()` or type-specific functions like `is.list()`.

```r
# Current deprecated pattern:
if(class(eqSS)!="list") stop("eqSS is not of class list")

# Recommended pattern:
if(!is.list(eqSS)) stop("eqSS is not of class list")
```

**Impact:** While not a performance issue per se, this pattern can cause issues with objects that have multiple classes and generates warnings in newer R versions.

## 3. Inefficient Loop-Based Range Generation

**Location:** `R/utils.R`, lines 61-69

**Function:** `ranges()`

**Issue:** Uses a for loop to create a list of ranges when this could be done more efficiently with `lapply` or vectorized operations.

```r
# Current pattern:
ranges <- function(diff, end){
  r <- list();
  for(i in 1:(end/diff)){
    r[[i]] <- (1 + ((i-1) * diff)) : (diff*i)
  }
  return(r)
}
```

**Recommended Fix:** Use `lapply` with `seq_len`:
```r
ranges <- function(diff, end){
  lapply(seq_len(end/diff), function(i) {
    (1 + ((i-1) * diff)) : (diff*i)
  })
}
```

## 4. Inefficient List Building in `flatten_barcode_vectors()`

**Location:** `R/stat_genetics_utils.R`, lines 317-343

**Issue:** Uses a for loop with incremental list index assignment. While not as bad as vector concatenation, pre-allocation and using `lapply` would be cleaner and potentially faster.

## 5. Repeated `which()` Calls in Matrix Adjustment

**Location:** `R/stat_genetics_utils.R`, lines 384-430

**Function:** `adjust_matrix_to_plaf()`

**Issue:** Uses a for loop over columns with repeated `which()` calls to find indices. For large matrices, this could be optimized using vectorized operations.

## 6. Nested Loops in COI Calculation

**Location:** `R/plotting.R`, lines 192-212

**Function:** `convert_barcode_vectors()`

**Issue:** Uses a for loop for COI calculation that could potentially be vectorized using `vapply` or `sapply`.

```r
# Current pattern:
for(i in 1:length(int.out)){
  if(n.strains[i]!=0){
    if(COI_type == "pcr_imperial"){
      COI[i] <- COI_pcr_imperial(i)
    }
    ...
  }
}
```

## Summary

The most impactful fix would be addressing issue #1 (the `ztrgeomintp` and `ztrnbinom` functions) as these are called during the main simulation pipeline and the inefficient vector growth pattern is a well-known R anti-pattern that can significantly slow down execution when many iterations are needed.
