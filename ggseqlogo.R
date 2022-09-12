install.packages("ggseqlogo")
# load required packages
require(ggplot2)
require(ggseqlogo)
# sample data
data(ggseqlogo_sample)
# plot logo
ggplot() + geom_logo( seqs_dna$MA0001.1 ) + theme_logo()
# shortcut
ggseqlogo( seqs_dna$MA0001.1 )
# plot based on probability
p2 = ggseqlogo( seqs_dna$MA0001.1, method = 'prob' )
p2
# can specify sequence type
ggseqlogo( seqs_aa$AKT1, seq_type='aa' )
# change color scheme to base-pairing
ggseqlogo(seqs_dna$MA0001.1, col_scheme='base_pairing')
# Create custom colour scheme
cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), groups=c('gr1', 'gr1', 'gr2', 'gr2'), 
                      cols=c('purple', 'purple', 'blue', 'blue'))

# Generate sequence logo
ggseqlogo(seqs_dna$MA0001.1, col_scheme=cs1)
# differnt scheme
cs2 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), values=1:4)
ggseqlogo(seqs_dna$MA0001.1, col_scheme=cs2)
# mulitple logos
ggseqlogo(seqs_dna, ncol=4)
# different fonts
fonts = list_fonts(F)

p_list = lapply(fonts, function(f){
  ggseqlogo(seqs_dna$MA0001.1, font=f) + ggtitle(f)
})

do.call(gridExtra::grid.arrange, c(p_list, ncol=2))
# annotate and label
ggplot() + 
  annotate('rect', xmin = 0.5, xmax = 3.5, ymin = -0.05, ymax = 1.9, alpha = .1, col='black', fill='yellow') +
  geom_logo(seqs_dna$MA0001.1, stack_width = 0.90) + 
  annotate('segment', x = 4, xend=8, y=1.2, yend=1.2, size=2) + 
  annotate('text', x=6, y=1.3, label='Text annotation') + 
  theme_logo()
## Combine plots
# Sequences we're going to use for the logo
seqs = seqs_dna$MA0008.1

# Generate the sequence logo
p1 = ggseqlogo(seqs) + theme(axis.text.x = element_blank())

# Make data for sequence alignment
aln = data.frame(
  letter=strsplit("AGATAAGATGATAAAAAGATAAGA", "")[[1]], 
  species = rep(c("a", "b", "c"), each=8),
  x       = rep(1:8, 3)
)
aln$mut = 'no'
aln$mut[ c(2,15,20,23) ] = 'yes'

# Generate the sequence alignment
p2 = ggplot(aln, aes(x, species)) +
  geom_text(aes(label=letter, color=mut, size=mut)) + 
  scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + xlab('') + 
  scale_color_manual(values=c('black', 'red')) + 
  scale_size_manual(values=c(5, 6)) + 
  theme_logo() + 
  theme(legend.position = 'none', axis.text.x = element_blank()) 

# Generate barplot data
bp_data = data.frame(x=1:8, conservation=sample(1:100, 8))

# Generate barplot data 
p3 = ggplot(bp_data, aes(x, conservation)) +
  geom_bar(stat='identity', fill='grey') + 
  theme_logo() + 
  scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + 
  xlab('')


# Now combine using cowplot, which ensures the plots are aligned
suppressMessages( require(cowplot) )
plot_grid(p1, p2, p3,  ncol = 1, align = 'v')

# Try with pcre data
# align mulitple cres
# read in spreadsheet with aligned_pcre\tmotif
nc= "M1235_1.01.fa.aln.mod.txt"

p2<-read.table(nc, header=T, sep='\t', row.names=1) 
p2l <- list(M1235_1.01=rownames(subset(p2, MOTIF == 'M1235_1.01')))
print(p2l)
# plot logo with probability score
ggseqlogo(p2l$M1235_1.01, method = 'prob',seq_type='DNA')
