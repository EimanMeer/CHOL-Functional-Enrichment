# --- PROJECT: FUNCTIONAL ENRICHMENT ANALYSIS ---
# PURPOSE: Identify Biological Pathways (GO/KEGG)
# AUTHOR: Eiman Meer

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 1. Prepare your gene list (using the sig_list from your previous script)
# We need Entrez IDs for clusterProfiler
gene_ids <- bitr(sig_list$gene_name, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# 2. Run Gene Ontology (GO) Enrichment
# This finds specific biological processes that are over-represented
ego <- enrichGO(gene          = gene_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP", # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                readable      = TRUE)

# 3. Visualization: The DotPlot
# This is a classic "Bioinformatics" figure for your GitHub
dotplot(ego, showCategory=20) + 
  ggtitle("Top Biological Processes in Cholangiocarcinoma")

# 1. Ensure the ego object is using Gene Symbols (makes the plot readable)
ego_readable <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

# 2. Create a STRICTLY matching fold-change vector
# We filter the sig_list to only include genes that actually made it into the enrichment
sig_matched <- sig_list %>% filter(gene_name %in% names(ego_readable@geneSets))

# Create the named numeric vector
genelist_for_plot <- sig_matched$log2FoldChange
names(genelist_for_plot) <- sig_matched$gene_name

# 3. Final Cnetplot with explicit size setting
# We set categorySize = "pvalue" to size the pathway dots by significance
cnetplot(ego_readable, 
         foldChange = genelist_for_plot, 
         showCategory = 5, 
         categorySize = "pvalue", # Sizing dots by p-value
         colorEdge = TRUE, 
         shadow_text = TRUE)
# 1. Ensure you have the latest mapping
ego_readable <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')

# 2. Reset the fold-change vector names to be sure
genelist_for_plot <- sig_list$log2FoldChange
names(genelist_for_plot) <- sig_list$gene_name

# 3. Simple Cnetplot (No extra arguments that cause crashes)
library(enrichplot)
cnetplot(ego_readable, 
         foldChange = genelist_for_plot, 
         showCategory = 5)
ggsave("CHOL_Gene_Concept_Network.png", plot = last_plot(), 
       width = 12, height = 8, dpi = 300)