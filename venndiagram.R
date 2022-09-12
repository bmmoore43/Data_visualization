library(VennDiagram)

#pairwise venn diagram example
pdf("GO-aracyc_PM_venn.pdf")
venn.plot <- draw.pairwise.venn(
  area1 = 1644, #gene number of A
  area2 = 2006, #gene number of B
  cross.area = 744, #overlap
  category = c("Aracyc", "GO"), #label
  fill = c("darkorchid1","green"), #fill
  lty = "blank",
  cex = 2, #text size & positions
  cat.cex = 1,
  cat.pos = c(275,110),
  cat.dist = c(0.01,0.01),
  cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 4,
  ext.dist = 1,
  ext.length = 0.15,
  ext.line.lwd = 2,
  ext.line.lty = "dashed");
grid.draw(venn.plot);
dev.off()

# A simple three-set diagram
venn.plot <- draw.triple.venn(65, 75, 85,
                              35, 15, 25, 5, c("First", "Second", "Third"));
grid.draw(venn.plot);
grid.newpage();

# A more complicated diagram
pdf("GO+aracyc_benchmark_enz_venn_PM.pdf")
venn.plot <- draw.triple.venn(
  area1 = 785,
  area2 = 45,
  area3 = 202,
  n12 = 2,
  n23 = 0,
  n13 = 0,
  n123 = 0,
  category = c("PM GO+Aracyc", "SM benchmark", "PM benchmark"),
  fill = c("darkorchid1", "goldenrod1", "green"),
  lty = "blank",
  cex = 2.5,
  cat.cex = 1.5,
  cat.dist = 0.008,
  cat.col = c("darkorchid4", "darkorange", "darkgreen")
);
grid.draw(venn.plot);
dev.off()

# quad venn diagram
pdf ("SMvsPM_tomatocyc-BM_venn.pdf")
venn.plot <- draw.quad.venn(
  area1 = 356, #tomatocyc SM
  area2 = 2261, #tomatocyc PM
  area3 = 68, #gold SM
  area4 = 13, #gold PM
  n12 = 85,
  n13 = 25,
  n23 = 25,
  n14 = 3,
  n24 = 9,
  n34 = 2,
  n123 = 8,
  n124 = 2,
  n134 = 1,
  n234 = 1,
  n1234 = 1,
  category = c("Tomatocyc SM", "Tomatocyc PM", "benchmark SM", "benchmark PM"),
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  lty = "blank",
  cex = 1,
  cat.cex = 1.1,
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.dist = 0.2
);

grid.draw(venn.plot);
dev.off()
