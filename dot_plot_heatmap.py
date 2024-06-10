import sys, os, re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime
from operator import itemgetter


def chunk_generator( seq, chunk_size):
    """! @brief generates successive n-sized chunks from input sequence """
    
    for i in xrange( 0, len( seq ), chunk_size):
        yield seq[ i:i + chunk_size]


def construct_seq_block_file( input_file, output_file, block_size ):
	"""! @brief construct file with chunked sequence """
	
	counter = 0
	with open( output_file, "w" ) as out:
		with open( input_file, "r" ) as f:
			f.readline() #remove header
			seq = []
			line = f.readline()
			while line:
				if line[0] == ">":
					chunks = chunk_generator( "".join( seq ), block_size)
					for each in chunks:
						out.write( '>' + str( counter ).zfill( 8 ) + '\n' + each + '\n' )
						counter += 1
					seq = []
				else:
					seq.append( line.strip() )
				line = f.readline()
			chunks = chunk_generator( "".join( seq ), block_size)
			for each in chunks:
				out.write( '>' + str( counter ).zfill( 8 ) + '\n' + each + '\n' )
				counter += 1


def load_blast_results( blast_result_file ):
	"""! @brief load blast results for plot """
	
	blast_results = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				hit = blast_results[ parts[0] ]
				if float( parts[-1] ) > hit['score']:
					del blast_results[ parts[0] ]
					blast_results.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } } )
			except KeyError:
				blast_results.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': ( int( parts[8] )+int( parts[9] ) )*0.5, 'score': float( parts[-1] ) } } )
			line = f.readline()
	return blast_results


def calculate_color( value ):
	"""! @brief calculates color for given value """
		
	r = 255-int( 255*value**4 )
	g = 255-int( 255*value**4 )
	b = 255
		
	color = '#%02x%02x%02x' % (r, g, b)
	
	return color


def construct_dot_plot( self_blast_results, other_blast_results, output_figure, show_status, namex, namey, title ):
	"""! @brief construct dot plot with score weighted points """
	
	# --- construction of plot --- #
	fig, ax = plt.subplots( figsize=(10,10) )
	
	maxx = 0
	maxy = 0
	for key in self_blast_results.keys():
		try:
			other_point = other_blast_results[ key ]
			self_point = self_blast_results[ key ]
			score = other_point[ 'score' ] / self_point[ 'score' ]
			point_color = calculate_color( score )
			x = self_point['pos'] / 1000.0
			if x > maxx:
				maxx = 0 + x
			y = other_point['pos'] / 1000.0
			if y > maxy:
				maxy = 0 + y
			ax.scatter( x, y, color=point_color,s=1 )
		except KeyError:
			pass
	
	ax.set_xlabel( namex + " [kbp]" )
	ax.set_ylabel( namey + " [kbp]" )
	
	ax.set_xlim( 0, maxx )
	ax.set_ylim( 0, maxy )
	
	ax.ticklabel_format( axis='y', scilimits=(-2,2) )	#, style="sci"
	ax.ticklabel_format( axis='x', scilimits=(-2,2) )	#, style="sci"
	
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	ax.get_yaxis().get_major_formatter().set_scientific(False)
	
	patches = []
	for i in range( 11 ):
		patches.append( mpatches.Patch(color=calculate_color( i/10.0 ), label=str(i/10.0)  ) )
	ax.legend( handles=patches, loc='upper left' )	#'lower right'
	
	ax.set_title( title )
	
	if show_status:
		plt.show()
	
	fig.savefig( output_figure, dpi=1200 )


def load_all_seqs_from_multiple_fasta_file( filename ):
	"""! @brief load all sequences from multiple fasta file """
	
	data = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				data.update( { header: "".join( seq ) } )
				header = line.strip()[1:].split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		data.update( { header: "".join( seq ) } )
	return data


def generate_concatenated_seq( infile, outfile ):
	"""! @brief generate new sequence file with concatenated sequence """
	
	seq = []
	with open( infile, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '>':
				seq.append( line.strip() )
			line = f.readline()
	
	with open( outfile, "w" ) as out:
		out.write( '>seq\n' + "".join( seq ) + '\n' )


def load_gene_positions( gff1_file ):
	"""! @brief load gene positions from given GFF file """
	
	gene_positions = {}
	with open( gff1_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					ID = re.findall( "AT\dG\d+", parts[-1] )[0]
					try:
						gene_positions[ parts[0] ].append( { 'start': int( parts[3] ), 'end': int( parts[4] ), 'agi': ID } )
					except KeyError:
						gene_positions.update( { parts[0]: [ { 'start': int( parts[3] ), 'end': int( parts[4] ), 'agi': ID } ] } )
			line = f.readline()
	return gene_positions


def load_BAC_positions( gff2_file ):
	"""! @brief load BAC IDs from GFF file """
	
	BAC_IDs = {}
	with open( gff2_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				try:
					BAC_IDs[ parts[0] ].append( { 'start': int( parts[3] ), 'end': int( parts[4] ), 'id': parts[-1] } )
				except KeyError:
					BAC_IDs.update( { parts[0]: [ { 'start': int( parts[3] ), 'end': int( parts[4] ), 'id': parts[-1] } ] } )
			line = f.readline()
	return BAC_IDs


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	in_seq_ref_file1 = arguments[  arguments.index( '--in1' )+1 ]
	in_seq_ref_file2 = arguments[  arguments.index( '--in2' )+1 ]
	
	if '--show' in arguments:
		show_status = True
	else:
		show_status = False
		
	if '--block' in arguments:
		block_size = int( arguments[  arguments.index( '--block' )+1 ] )
	else:
		block_size = 1000
	
	if '--namex' in arguments:
		namex = arguments[  arguments.index( '--namex' )+1 ]
	else:
		namex = "x"
	
	if '--namey' in arguments:
		namey = arguments[  arguments.index( '--namey' )+1 ]
	else:
		namey = "y"
	
	if '--name' in arguments:
		name = arguments[  arguments.index( '--name' )+1 ]
	else:
		name = ""

	if '--title' in arguments:
		title = arguments[  arguments.index( '--title' )+1 ]
	else:
		title = ""
	
	if '--gff1' in arguments:
		gff1_file = arguments[  arguments.index( '--gff1' )+1 ]
		genes = load_gene_positions( gff1_file )
	else:
		genes = {}
	
	if '--gff2' in arguments:
		gff2_file = arguments[  arguments.index( '--gff2' )+1 ]
		BACs = load_BAC_positions( gff2_file )
	else:
		BACs = {}
	
	if '--chr' in arguments:
		chromosome = arguments[  arguments.index( '--chr' )+1 ]
	else:
		chromosome = ""
	
	if '--start' in arguments:
		start = int( arguments[  arguments.index( '--start' )+1 ] )
	else:
		start = 0
	
	if '--end' in arguments:
		end = int( arguments[  arguments.index( '--end' )+1 ] )
	else:
		end = 0
	
	# --- construct complex name --- #
	try:
		genes_per_chr = sorted( genes[ chromosome ], key=itemgetter('start') )
		gene_names = []
		for gene in genes_per_chr:
			if gene['start'] < end:
				if gene['end'] > start:
					gene_names.append( gene['agi'] )
		if len( gene_names ) > 0:
			if len( gene_names ) > 2:
				title += "_" + gene_names[0] + "..." + gene_names[-1]
			else:
				title += "_" + ";".join( gene_names )
	except KeyError:
		pass
	
	try:
		bacs_per_chr = sorted( BACs[ chromosome ], key=itemgetter('start') )
		bac_names = []
		for bac in bacs_per_chr:
			if bac['start'] < end:
				if bac['end'] > start:
					bac_names.append( bac['id'] )
		if len( bac_names ) > 0:
			if len( bac_names ) > 2:
				title += "_" + bac_names[0] + "..." + bac_names[-1]
			else:
				title += "_" + ";".join( bac_names )
	except KeyError:
		pass
	
	
	prefix = arguments[  arguments.index( '--out' )+1 ]
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	output_figure = prefix + name +  "_dot_plot_heatmap.pdf"
	seq_ref_file1 = prefix + name + "_seq1.fasta"
	seq_ref_file2 = prefix + name +  "_seq2.fasta"
	
	generate_concatenated_seq( in_seq_ref_file1, seq_ref_file1 )
	generate_concatenated_seq( in_seq_ref_file2, seq_ref_file2 )
	
	
	active = True
	
	# --- set output files --- #
	self_blast_result_file = prefix + "self_blast_result_file.txt"
	other_blast_result_file = prefix + "other_blast_result_file.txt"
	
	if active:
		seq_block_file1 = prefix + "seq_blocks.fasta"
		construct_seq_block_file( seq_ref_file1, seq_block_file1, block_size )
		
		# --- run blast vs. self and vs. other --- #
		self_blast_db = prefix + "self_blast_db"
		other_blast_db = prefix + "other_blast_db"
		if active:
			os.popen( "makeblastdb -in " + seq_ref_file1 + " -out " +  self_blast_db + " -dbtype nucl" )
			os.popen( "makeblastdb -in " + seq_ref_file2 + " -out " +  other_blast_db + " -dbtype nucl" )
		
		if active:
			os.popen( "blastn -query " + seq_block_file1 + " -db " + self_blast_db + " -out " +  self_blast_result_file + " -outfmt 6 -evalue 0.01 -num_threads 8" )
			os.popen( "blastn -query " + seq_block_file1 + " -db " + other_blast_db + " -out " +  other_blast_result_file + " -outfmt 6 -evalue 0.01 -num_threads 8" )
		
	# --- load blast results --- #
	self_blast_results = load_blast_results( self_blast_result_file )
	other_blast_results = load_blast_results( other_blast_result_file )
	
	
	print "Generating dot plot heatmap ..."
	construct_dot_plot( self_blast_results, other_blast_results, output_figure, show_status, namex, namey, title )


if '--cite' in sys.argv:
	sys.exit( __reference__ )
	
if '--help' in sys.argv or '-h' in sys.argv:
	sys.exit( __usage__ )

elif '--out' in sys.argv and '--in1' in sys.argv and '--in2' in sys.argv:
	main( sys.argv )

else:
	sys.exit( __usage__ )
