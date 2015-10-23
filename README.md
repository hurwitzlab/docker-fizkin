# Docker Fizink

Build or grab Jellyfish binaries:

    http://www.genome.umd.edu/jellyfish.html#Release

Install needed Perl modules here:

    $ cpanm -l . --self-contained --installdeps .

Build Docker image:

    $ docker build -t fizkin .

Run with args:

    $ docker run --rm -v /path/to/data:/work -w /work fizkin \
      run-fizkin --input_dir /work/fasta --output_dir /work/out

# See Also

    http://github.com/hurwitzlab/fizkin

# Authors

* Bonnie Hurwitz <bhurwitz@email.arizona.edu>
* Ken Youens-Clark <kyclark@email.arizona.edu>
