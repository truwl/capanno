#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [gawk]
stdout: $(inputs.outFileName)
requirements:
- class: InlineJavascriptRequirement

inputs:

  fieldSeparator:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -F
    doc: |
      --field-separator fs
              Use fs for the input field separator (the value of the FS predefined variable).

  assign:
    type: ["null", string]
    inputBinding:
      position: 1
      separate: false
      prefix: -v var=
    doc: |
      --assign var=val
              Assign the value val to the variable var, before execution of the program begins.  Such variable values are available to the BEGIN rule of an AWK program.

  charactersAsBytes:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -b
    doc: |
      --characters-as-bytes
              Treat all input data as single-byte characters. In other words, don't pay any attention to the locale information when attempting to process strings as  multibyte  characters.   The
              --posix option overrides this one.

  traditional:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -c
    doc: |
      --traditional
              Run in compatibility mode.  In compatibility mode, gawk behaves identically to Brian Kernighan's awk; none of the GNU-specific extensions are recognized.  See GNU EXTENSIONS, below for more information.

  copyright:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -C
    doc: |
      --copyright
              Print the short version of the GNU copyright information message on the standard output and exit successfully.

  dumpVariables:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -d
    doc: |
      -d[file]
       --dump-variables[=file]
              Print a sorted list of global variables, their types and final values to file.  If no file is provided, gawk uses a file named awkvars.out in the current directory.
              Having a list of all the global variables is a good way to look for typographical errors in your programs.  You would also use this option if you have a large program with a lot  of
              functions, and you want to be sure that your functions don't inadvertently use global variables that you meant to be local.  (This is a particularly easy mistake to make with simple
              variable names like i, j, and so on.)

  debug:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -D
    doc: |
      -D[file]
       --debug[=file]
              Enable debugging of AWK programs.  By default, the debugger reads commands interactively from the keyboard (standard input).  The optional file argument specifies a file with a list
              of commands for the debugger to execute non-interactively.

  programText:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -e
    doc: |
      -e program-text
       --source program-text
              Use  program-text as AWK program source code.  This option allows the easy intermixing of library functions (used via the -f and --file options) with source code entered on the command line.  It is intended primarily for medium to large AWK programs used in shell scripts.

  exec:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: -E
    doc: |
      -E file
       --exec file
              Similar to -f, however, this is option is the last one processed.  This should be used with #!  scripts, particularly for CGI applications, to avoid passing  in  options  or  source code (!) on the command line from a URL.  This option disables command-line variable assignments.

  genPot:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -g
    doc: |
      -g
       --gen-pot
              Scan  and parse the AWK program, and generate a GNU .pot (Portable Object Template) format file on standard output with entries for all localizable strings in the program.  The program itself is not executed.  See the GNU gettext distribution for more information on .pot files.

  help:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -h
    doc: |
      -h
       --help Print a relatively short summary of the available options on the standard output.  (Per the GNU Coding Standards, these options cause an immediate, successful exit.)

  includeFile:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -i
    doc: |
       -i include-file
       --include include-file
              Load an awk source library.  This searches for the library using the AWKPATH environment variable.  If the initial search fails, another attempt will be  made  after  appending  the .awk suffix.  The file will be loaded only once (i.e., duplicates are eliminated), and the code does not constitute the main program source.

  load:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      -l lib
       --load lib
              Load  a  shared library lib.  This searches for the library using the AWKLIBPATH environment variable.  If the initial search fails, another attempt will be made after appending the default shared library suffix for the platform.  The library initialization routine is expected to be named dl_load().

  lint:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -L
    doc: |
      -L [value]
       --lint[=value]
              Provide warnings about constructs that are dubious or non-portable to other AWK implementations.  With an optional argument of fatal, lint warnings become fatal errors.  This may be drastic,  but  its  use will certainly encourage the development of cleaner AWK programs.  With an optional argument of invalid, only warnings about things that are actually invalid are issued. (This is not fully implemented yet.)

  bigNum:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -M
    doc: |
      -M
       --bignum
              Force arbitrary precision arithmetic on numbers. This option has no effect if gawk is not compiled to use the GNU MPFR and MP libraries.

  nonDecimalData:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -n
    doc: |
      -n
       --non-decimal-data
              Recognize octal and hexadecimal values in input data.  Use this option with great caution!

  UseLcNumeric:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -N
    doc: |
      -N
       --use-lc-numeric
              This forces gawk to use the locale's decimal point character when parsing input data.  Although the POSIX standard requires this behavior, and  gawk  does  so when  --posix  is  in effect,  the default is to follow traditional behavior and use a period as the decimal point, even in locales where the period is not the decimal point character.  This option overrides the default behavior, without the full draconian strictness of the --posix option.

  prettyPrint:
    type: ["null", string]
    inputBinding:
      position: 1
      separate: false
      prefix: -o
    doc: |
       -o[file]
       --pretty-print[=file]
              Output a pretty printed version of the program to file.  If no file is provided, gawk uses a file named awkprof.out in the current directory.

  optimize:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -O
    doc: |
      -O
       --optimize
              Enable optimizations upon the internal representation of the program.  Currently, this includes simple constant-folding, and tail call elimination for recursive functions. The  gawk maintainer hopes to add additional optimizations over time.

  posix:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -P
    doc: |
      -P
       --posix
              This turns on compatibility mode, with the following additional restrictions:

              · \x escape sequences are not recognized.

              · Only space and tab act as field separators when FS is set to a single space, newline does not.

              · You cannot continue lines after ?  and :.

              · The synonym func for the keyword function is not recognized.

              · The operators ** and **= cannot be used in place of ^ and ^=.

  reInterval:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      -r
       --re-interval
              Enable the use of interval expressions in regular expression matching (see Regular Expressions, below).  Interval expressions were not traditionally available in the  AWK  language.
              The POSIX standard added them, to make awk and egrep consistent with each other.  They are enabled by default, but this option remains for use with --traditional.

  sandbox:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -S
    doc: |
      -S
       --sandbox
              Runs gawk in sandbox mode, disabling the system() function, input redirection with getline, output redirection with print and printf, and loading dynamic extensions.  Command execution (through pipelines) is also disabled.  This effectively blocks a script from accessing local resources (except for the files specified on the command line).

  lintOld:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -t
    doc: |
       -t
       --lint-old
              Provide warnings about constructs that are not portable to the original version of UNIX awk.

  version:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -V
    doc: |
      -V
      --version
              Print version information for this particular copy of gawk on the standard output.  This is useful mainly for knowing if the current copy of gawk on your system is up to  date  with respect  to  whatever  the Free Software Foundation is distributing.  This is also useful when reporting bugs.  (Per the GNU Coding Standards, these options cause an immediate, successful exit.)

  end:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --
    doc: |
       --     Signal the end of options. This is useful to allow further arguments to the AWK program itself to start with a “-”.  This provides consistency with the argument  parsing  convention used by most other POSIX programs.


  awkScript:
    type: ["null", File]
    inputBinding:
      position: 3
      prefix: -f
    doc: >
      Read the AWK program source from this program-file, instead of from the first command line argument.
      Multiple -f (or --file) options may be used.

  rawPatternAction:
    type: ["null", string]
    inputBinding:
      position: 4
    doc: String to specify awk instructions.

  inputFiles:
    type:
      - "null"
      - type: array
        items: File
    inputBinding:
      position: 5
    doc: Input files supplied to the awkFile.

  outFileName:
    type: string
    default: gawkOutput.txt
    doc: specify the name of the file output by gawk. This is read by $(stdout)


outputs:
  awkOutput:
    type: stdout
    doc: The output of the awk script.
