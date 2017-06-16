## Genome Assembly
	
  The	Assembly	Algorithm	implemented	in	Proj_1.py	makes	use	of	a	Euclidean	graph,	
described	by	a	nested	hash	table,	to	assemble	30	base	reads.	Upon	startup,	reads	are	split	into	
15	base	substrings,	the	first	15	reads	constituting	a	prefix,	and	the	last		15	a	suffix.	As	reads	are	
input,	line	by	line	from	a	.fasta	file,	both	prefixes	and	suffixes	are	encoded	to	a	uniquely	
identifying	index	number	where	each	nucleotide	(A	C	T	G),	is	denoted	by	an	integer,	zero,	one,	
two,	or	three,		respectively,	multiplied	by	four	raised	to	its	corresponding	position	in	the	
prefix/suffix	then	added	to	the	value	of	its	neighbor.	Positions	are	numbered	from	zero	to	
fourteen,	thus,	as	an	example,	the	string	ACTGACTGACTGACT	will	be	represented	by	the	index:	
(0)4 0 +(1)4 1 +(2)4 2 +(3)4 3 +(0)4 4 (1)4 5 +(2)4 6 +(3)4 7 +(0)4 8 +(1)4 9 +(2)4 10 (3)4 11 +(0)4 12 +(1)4 13 +(2)4 14 .	
Once	the	prefix	and	suffix	obtained	from	a	read	are	encoded	to	an	index	value,	the	prefix	is	
compared	to	the	stored	values	in	the	hash	table.	If	the	prefix	has	been	seen	before,	its	suffix	is	
stored	within	a	hash	table	of	suffixes	nested	within	a	hash	table	of	prefixes,	where	each	prefix	
acts	as	a	key	to	the	inner	hash	table	of	corresponding	suffixes.	The	dictionary	containing	the	
suffix	as	a	key	is	then	given	a	value	corresponding	to	the	count	of	the	reads	that	have	been	
read.	This	value	matches	the	key	to	another	dictionary	which	stores	the	count	as	key	with	the	
alphabetic	representation	of	the	read	as	a	value.		If	both		prefix	and	suffix	have	been	observed	
in	previous	iterations	of	a	.fasta	files	parse,	the	read	count	stored	in	the	suffix	of	the	inner	hash	
is	concatenated	with	the	read	number	of	the	previous	observation.	If	the	prefix,	suffix	pair	has	
not	been	observed,	a	new	prefix	is	added	to	the	outer	hash	along	with	a	new	suffix	to	the	inner	
hash,	referring	to	a	new	read	count	value.		As	reads	are	hashed	to	their	appropriate	location	in	
the	nested	hash	table,	each	prefix/suffix	is	given	a	out	degree	and	an	in	degree	where	the	out	
degree	refers	to	the	number	of	suffixes	that	issue	from	the	outer	prefix	of	a	hashed	value,	and	
the	in	degree	refers	to	the	number	of	times	that	a	hashed	value	from	the	outer	prefix	table	is	
represented	as	suffix	in	the	inner	hash	table.	These	steps	in	assembly	are	accomplished	by	the	
function	load_reads	in	Proj_1.py.	The	load_reads	function	returns	four	values.	A	nested	hash	
table	called,	hash_table	containing	the	outer	prefix	hash,	inner	prefix	hash	and	corresponding	
edge	numbers	represented	by	the	read	count,	two	hash	tables	read_table	and	
alpha_read_table,	containing	the	read	count	as	key	and	corresponding	prefix/suffix	value	as	the	
value,	and	degree_table,	a	hash	table	containing	a	prefix	value	as	key,	along	with	a	two	element	
list	containing	the	out	degree	and	in	degree	as	a	value.

The	values	returned	by	the	load	reads	function	are	used	to	initialize	attributes	within	a	
class	object	called	Path_Assembler.	The	purpose	of	using	a	class	to	store	the	hash	tables	
returned	by	the	load	reads	function	is	to	avoid	the	computational	cost	that	may	be	incurred	by	
repetitively	passing	the	sizable	nested	hash	table	that	results	from	load_reads	as	the	parameter	
of	whatever	function	may	be	used	to	find	paths	within	the	graph.	Python	uses	pass	by	
reference	value	to	pass	parameters	to	functions.	This	requires	making	a	copy	of	a	parameter	
each	time	it	is	passed	into	a	function.	This	drawback	in	using	python	can	be	avoided	by		housing	
the	nested	hash	table	data	structure	in	an	object	with	methods	to	search	the	hash	table.	This	
allows	the	graph	represented	by	the	nested	hash	table	to	be	searched	without	making	a	copy	every	time	the	graph	is	needed.	The	Hash	table	only	needs	to	be	copied	into	main	memory	
once,	using	a	class.	In	Proj_1.py	the	class	attribute	GraphHash	stores	the	nested	hash	table,	
ReadHash	stores	the	hash	table	containing	the	read	count	as	a	key	and	read	codes	as	values,	
DegreeHash	contains	the	hash	table	acting	as	a	reference	for	in	degrees	and	out	degrees,	and	
AlphaReadHash	contains	the	hash	table	mapping	read	counts	to	the	alphabetic	representation	
of	the	prefixes	and	suffixes.

Once	the	reads	are	loaded	into	a	nested	hash	table	and	stored	within	a	class,	the	
resulting	hash	table	may	be	searched	as	an	adjacency	list,	representing	a	graph,	to	find	
connecting	paths	representing	sections	of	a	genome.	Successive	searches	are	made	by	starting	
at	a	prefix	with	either	an	in	degree	greater	than	out	degree,	or,	if	no	such	prefix	is	available,	an	
in	degree	equal	to	out	degree,	where	the	out	degree	is	not	equal	to	zero.	In	the	Proj_1.py	
program,	the	method	start_node,	handles	the	task	of	finding	a	suitable	candidate	in	the	
GraphHash	nested	hash	table.	It	returns	either	the	index	value	of	a	proper	start	node,	or	-1,	for	
the	case	that	no	such	start	node	exists.

The	method	start_assembly,	iteratively	calls	start_node,	and	passes	the	key	generated	
by	the	start_node	method	to	the	method	find	paths.	Find	paths	is	responsible	for	traversing	the	
graph,	as	well	as	finding	and	storing	paths	within	the	class	attribute	paths.	

Once	a	start	node	is	passed	to	the	find_paths	method,	the	find_paths	method	begins	by	
checking	the	number	of	values	stored	within	the	start	node,	if	more	than	zero	suffixes	are	
stored	within	the	method,	the	suffixes	contain	edges,	and	the	outdegree	of	the	suffix	is	greater	
than	zero,	the	method	proceeds.	After	satisfying	the	initial	conditions,	find_paths	chooses	a	
suffix	stored	within	the	inner	hash	table	of	the	start	node	prefix.	This	is	accomplished	by	the	
method	suffix_picker.	Suffix	picker	chooses	a	suffix	that	has	not	been	used.	If	the	case	arises	
where	the	suffix	is	equal	to	the	prefix,	a	cycle	has	been	found,	the	in	degree/out	degree	is	
updated,		and	the	resulting	suffix/prefix	is	added	to	the	paths	as	a	lone	path.	Cycles	are	always	
avoided	by	the	conditions	imposed	by	suffix	picker,	and	added	to	the	paths	attribute	only	when	
no	other	options	are	available.	The	value	returned	by	suffix	picker	is	stored	in	rindex	within	
find_paths.	Once	found	the	edge	is	taken	out	of	the	GraphHash	data	structure,	indegree	and	
outdegrees	are	updated,	and	the	start	node	and	start	node	suffix	pass	through	a	set	of	
conditions	to	check	for	cycles	or	the	end	of	a	path.	If	either	the	end	of	the	path	or	an	impending	
cycle	is	detected,	the	current	path	traveled	is	stored	in	the	paths	attribute,	and	the	
start_assembly	method	makes	a	call	to	start_node	to	start	traversing	a	new	path.	If	neither	of	
these	conditions	arises,	the	find_paths	method	receives	a	recursive	call	to	begin	again	with	the	
suffix	as	the	new	startnode.	The	methods	start_assembly,	start_node,	and	find_paths	are	called	
until	there	are	no	longer	any	edges	to	traverse	in	the	graph.	

After	a	path	is	exhausted,	the	edges	that	have	been	used	to	traverse	the	graph	are	
stored	within	a	list	composed	of	lists	representing	paths.	The	searched	paths	are	removed	from	
the	nested	hash	table,	and	another	search	is	initiated.	This	process	continues	until	all	edges	and	
nodes	have	been	removed	from	the	graph,	and	there	is	nothing	left	to	search.	What	results	is	a	
list	of	paths	that	may	be	interpreted	as	contigs.	These	contigs	can	be	assembled	to	represent	
the	genome	of	whatever	organism’s	genetic	data	has	been	passed	to	the	Proj_1.py	program.	

To	assemble	the	contigs,	first,	all	duplicates	are	removed	by	applying	the	function	
remove_dup	to	the	paths	attribute.	Next	each	path	within	the	paths	attribute	is	merged	to	
create	a	sequence	of	characters	representing	each	path	with	the	combine_paths	function.	The	combine_paths	function	also	sorts	the	contigs	in	order	of	decreasing	length.	The	result	is	stored	
in	the	contigs	attribute.	After	the	lists	containing	paths	has	been	combined,	any	contig	that	is	
contained	within	another	contig	is	removed	with	the	remove_contained_contig	method,	
updating	the	contig	attribute	to	contain	only	contigs	that	cannot	be	placed	entirely	within	other	
contigs.	For	yeast	this	reduces	the	number	of	contigs	by	half,	from	around	2000	to	aproximately	
1000.	The	next	step	in	the	process	is to	compare	all	the	remaining	contigs,	find	
overlaps,	and	combine	the	results.
