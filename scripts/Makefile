sarah:
	LD_LIBRARY_PATH=/usr/local/imicrobe/docker-fizkin/local/lib ./run-fizkin.pl -i ~/work/sarah/fasta -o ~/work/sarah/fizkin --max_samples=300000 --debug

virome_t2: 
	LD_LIBRARY_PATH=/usr/local/imicrobe/docker-fizkin/local/lib ./run-fizkin.pl -i ~/work/skin/fasta/male28_virome -o ~/work/skin/male28_virome_t2 --metadata ~/work/skin/fasta/male28_virome.meta --max_samples=300000 --files=Ac2.fa,Ax2.fa,Fh2.fa,Pa2.fa,Ra2.fa,Sc2.fa,Tw2.fa,Um2.fa

virome_t3: 
	LD_LIBRARY_PATH=/usr/local/imicrobe/docker-fizkin/local/lib ./run-fizkin.pl -i ~/work/skin/fasta/male28_virome -o ~/work/skin/male28_virome_t3 --metadata ~/work/skin/fasta/male28_virome.meta --max_samples=300000 --files=Ac3.fa,Ax3.fa,Fh3.fa,Pa3.fa,Ra3.fa,Sc3.fa,Tw3.fa,Um3.fa

microbiome_t2: 
	LD_LIBRARY_PATH=/usr/local/imicrobe/docker-fizkin/local/lib ./run-fizkin.pl -i ~/work/skin/fasta/male28_microbiome -o ~/work/skin/male28_microbiome_t2 --metadata ~/work/skin/fasta/male28_microbiome.meta --max_samples=100000 --files Ac2.fa,Ax2.fa,Fh2.fa,Pa2.fa,Ra2.fa,Sc2.fa,Tw2.fa,Um2.fa

microbiome_t3: 
	LD_LIBRARY_PATH=/usr/local/imicrobe/docker-fizkin/local/lib ./run-fizkin.pl -i ~/work/skin/fasta/male28_microbiome -o ~/work/skin/male28_microbiome_t3 --metadata ~/work/skin/fasta/male28_microbiome.meta --max_samples=100000 --files Ac3.fa,Ax3.fa,Fh3.fa,Pa3.fa,Ra3.fa,Sc3.fa,Tw3.fa,Um3.fa

male_viromes: 
	LD_LIBRARY_PATH=/usr/local/imicrobe/docker-fizkin/local/lib ./run-fizkin.pl -i ~/work/skin/fasta/male_viromes -o ~/work/skin/male_viromes --metadata ~/work/skin/fasta/male_viromes.meta --max_samples=100000
