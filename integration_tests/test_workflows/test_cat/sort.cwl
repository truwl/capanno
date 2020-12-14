
cwlVersion: v1.0
class: CommandLineTool
baseCommand: sort
stdout: $(inputs.outputFile)
# stdin: need to put in some javascript logic here. If an input file is not specified, read stdin/
label: |
        sort [OPTION]... [FILE]...
        sort [OPTION]... --files0-from=F

inputs:

  ignoreLeadingBlanks:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -b
    doc: |
       -b, --ignore-leading-blanks
              ignore leading blanks


  dictionaryOrder:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -d
    doc: |
      -d, --dictionary-order
              consider only blanks and alphanumeric characters

  ignoreCase:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -f
    doc: |
       -f, --ignore-case
              fold lower case to upper case characters

  generalNumericSort:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -g
    doc: |
       -g, --general-numeric-sort
              compare according to general numerical value

  ignoreNonprinting:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -i
    doc: |
      -i, --ignore-nonprinting
              consider only printable characters

  monthSort:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -M
    doc: |
      -M, --month-sort
              compare (unknown) < 'JAN' < ... < 'DEC'

  humanNumericSort:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -h
    doc: |
         -h, --human-numeric-sort
              compare human readable numbers (e.g., 2K 1G)

  numericSort:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -n
    doc: |
       -n, --numeric-sort
              compare according to string numerical value
              compare according to string numerical value-S

  randomSort:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -R
    doc: |
       -R, --random-sort
              shuffle, but group identical keys.  See shuf(1)

  randomSource:
    type: ["null", File]
    inputBinding:
      position: 1
      separate: false
      prefix: --random-source=
    doc: |
      --random-source=FILE
              get random bytes from FILE

  reverse:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -r
    doc: |
        -r, --reverse
              reverse the result of comparisons

  wordSort:
    type: ["null", string]
    inputBinding:
      position: 1
      separate: false
      prefix: --sort=
    doc: |
      --sort=WORD
              sort according to WORD: general-numeric -g, human-numeric -h, month -M, numeric -n, random -R, version -V

  versionSort:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -V
    doc: |
      -V, --version-sort
              natural sort of (version) numbers within text

# Other options

  batchSize:
    type: ["null", int]
    inputBinding:
      position: 2
      separate: false
      prefix: --batch-size=
    doc: |
       --batch-size=NMERGE
              merge at most NMERGE inputs at once; for more use temp files

  check:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -c
    doc: |
      -c, --check, --check=diagnose-first
              check for sorted input; do not sort

  quietCheck:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: -C
    doc: |
      -C, --check=quiet, --check=silent
              like -c, but do not report first bad line

  compressProgram:
    type: ["null", string]
    inputBinding:
      position: 1
      separate: false
      prefix: --compress-program=
    doc: |
      --compress-program=PROG
              compress temporaries with PROG; decompress them with PROG -d

  debug:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --debug
    doc: |
      --debug
              annotate the part of the line used to sort, and warn about questionable usage to stderr

  key:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -k
    doc: |
      -k, --key=KEYDEF
              sort via a key; KEYDEF gives location and type

  merge:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      -m, --merge
              merge already sorted files; do not sort

  outputFile:
    type: string
    default: sortOutput.txt
#    inputBinding:
#      position: 2
#      prefix: -o
#    doc: |
#      -o, --output=FILE
#              write result to FILE instead of standard output

  stable:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: -s
    doc: |
      -s, --stable
              stabilize sort by disabling last-resort comparison

  bufferSize:
    type: ["null", string]
    inputBinding:
      position: 2
      prefix: -S
    doc: |
       -S, --buffer-size=SIZE
              use SIZE for main memory buffer. e.g. 60G

  fieldSeparator:
    type: string?
    inputBinding:
      position: 2
      prefix: -t
    doc: |
      -t, --field-separator=SEP
              use SEP instead of non-blank to blank transition

  temporaryDirectory:
    type: ["null", Directory]
    inputBinding:
      position: 2
      prefix: -T
    doc: |
       -T, --temporary-directory=DIR
              use DIR for temporaries, not $TMPDIR or /tmp; multiple options specify multiple directories

  parallel:
    type:  ["null", int]
    inputBinding:
      position: 2
      separate: false
      prefix: --parallel=
    doc: |
      --parallel=N
              change the number of sorts run concurrently to N

  unique:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: -u
    doc: |
      -u, --unique
              with -c, check for strict ordering; without -c, output only the first of an equal run

  zeroTerminated:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: -z
    doc: |
     -z, --zero-terminated
              line delimiter is NUL, not newline

  help:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --help
    doc: |
      --help display this help and exit

  version:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --version
    doc: |
      --version
              output version information and exit

# specify input

  files0From:
    type:  ["null", File]
    inputBinding:
      position: 3
      separate: false
      prefix: --files0-from=
    doc: |
      --files0-from=F
              read input from the files specified by NUL-terminated names in file F; If F is - then read names from standard input

  inputFile:
    type: ["null", File]
    inputBinding:
      position: 3
    doc: Specify files to sort.

#  stdInput:
#    type: string?
#    inputBinding:
#      position: 3
#    doc: Wnen no files are specified, read standard input.



outputs:
  output:
    type: stdout
