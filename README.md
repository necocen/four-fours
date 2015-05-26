# four-fours' solver

## Build
* requires C++11
* requires boost(boost/optional, boost/program_options, `-lboost_program_options`)
* (OPTIONAL) requires OpenMP(`-fopenmp`) for concurrent computing

## Usage

```sh
# solve normal 4,4,4,4 four-fours and print results on 4444.html
$ four_fours
# solve 1,2,3,4
$ four_fours 1234
```

```
Options:
  -h [ --help ]                 Shows this help
  -u [ --unary-limit ] arg (=2) Maximum number of consecutive application of
                                unary operators.
  -q [ --quiet ] arg            Suppress progress reporting.
  --no-cutoff-sqrt              Do not restrict sqrt operand.
  -o [ --output-file ] arg      Output file name.
```

## Acknowledgement
* Using Boost.Optional, Boost.Program_options
* Using MathJax for visualizing equations
