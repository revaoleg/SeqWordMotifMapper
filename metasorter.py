import sys, os
import main, tools, cui, seq_io

########################################################################
class MetaSorter:
    def __init__(self, version="", date_of_creation=""):
        self.version = version
        self.date_of_creation = date_of_creation
    
        # Input/output
        self.oIO = seq_io.IO()
        

    def parse_gff_file(self, gff_path, contigs):
        # Parse GFF file
        if not gff_path or not os.path.exists(gff_path):
            raise ValueError(f"Input GFF file {gff_path} does not exists!")
        tools.msg(f"Reading {gff_path}")
        dGFF = self.oIO.readGFF(path=gff_path, contigs=contigs, mode='dictionary')
        dGFF['input file'] = gff_path
        return [dGFF]
        
    def select_bin_gff(self, master_gff, contigs):
        records = [rec for rec in master_gff if rec.startswith(tuple(contigs))]
        return records
    
    def group_contigs(self, contig_titles, group_collection, unresolved_contigs="keep"):
        # If unresolved-contigs == 'keep', unmethylated contigs kept binned, if 'split', separated to individual bins 
        if len(group_collection) == 0 and unresolved_contigs == "keep":
            return [contig_titles]
            
        def split_overlapping_groups(raw_groups):
            """Split partly overlapping lists into non-overlapping consecutive groups."""
            title_to_mask = {}
    
            for i, group in enumerate(raw_groups):
                for title in group:
                    title_to_mask.setdefault(title, set()).add(i)
    
            result = []
            current_group = []
            current_mask = None
    
            for title in contig_titles:
                if title not in title_to_mask:
                    continue
    
                mask = frozenset(title_to_mask[title])
    
                if current_mask is None or mask != current_mask:
                    if current_group:
                        result.append(current_group)
                    current_group = [title]
                    current_mask = mask
                else:
                    current_group.append(title)
    
            if current_group:
                result.append(current_group)
    
            return result
    
        groups = []
        used = set()
    
        # 1. Positive groups
        positive_raw = [
            record.get("positive", [])
            for record in group_collection
            if record.get("positive")
        ]
    
        positive_groups = split_overlapping_groups(positive_raw)
        groups.extend(positive_groups)
    
        for group in positive_groups:
            used.update(group)
    
        # 2. Negative groups, excluding positive titles
        negative_raw = [
            [title for title in record.get("negative", []) if title not in used]
            for record in group_collection
            if record.get("negative")
        ]
    
        negative_groups = split_overlapping_groups([g for g in negative_raw if g])
        groups.extend(negative_groups)
    
        for group in negative_groups:
            used.update(group)
    
        # 3. Median titles, excluding positive and negative titles
        median_titles = []
    
        for record in group_collection:
            for title in record.get("median", []):
                if title not in used and title not in median_titles:
                    median_titles.append(title)
    
        if median_titles:
            median_group = [title for title in contig_titles if title in median_titles]
            groups.append(median_group)
            used.update(median_group)
    
        # 4. Titles absent from all grouped categories
        for title in contig_titles:
            if title not in used:
                groups.append([title])
                used.add(title)
    
        return groups
         
    def execute(self, args):       
        # Set runtime options
        options = cui.get_options()
        options['-u'] = args.input_folder               # Input folder
        options['-o'] = args.output_folder              # Output folder
        options['-d'] = args.project_folder             # Project folder
        options['-d'] = args.project_folder             # Project folder
        options['-c'] = args.mqs                        # Minimum methylation/modification quality score
        options['-ogf'] = args.output_graph_format      # Graphical format
        options['-dpi'] = args.dpi                      # DPI for rastr format
        options['-wl'] = args.sliding_window_length     # Sliding window length
        options['-ws'] = args.sliding_window_step       # Sliding window step
        options['-mm'] = args.circular_graph.upper()    # Request circular graph
        options['-sp'] = "Y"                            # Request statistic panel
        options['-s'] = "sites"                         # Request site statistics instead of motif statistics
        options['-f'] = "M"                             # Request statistics of methylated sites instead of unmethylated 'U'
        options['-r'] = 0                               # No limits for number of validated sites (server application)
        
        separator = args.file_name_separator
        
        # List of motifs
        motifs = [s.strip() for s in args.motifs.split(";")]
        
        if len(motifs) == 0:
            sys.exit("No motifs were specified!")
            
        fasta_extensions = (".FA", "FASTA", ".FNA", ".FST")
        genbank_extensions = (".GB", ".GBK", ".GBF", ".GBFF")
        # Current version works only with FASTA files
        genbank_extensions = ()
        
        input_path = os.path.join(args.input_folder, args.project_folder)
        master_gff = None
        if args.gff_file:
            master_gff_file_path = os.path.join(input_path, args.gff_file)
            if os.path.exists(master_gff_file_path):
                master_gff = self.oIO.open_text_file(path = master_gff_file_path, flg_inlist = True)
        
        # Loop through pairs of GFF and sequence files
        for seq_fname in [fn for fn in os.listdir(os.path.abspath(os.path.join(args.input_folder, args.project_folder))) 
                if fn.upper().endswith(fasta_extensions + genbank_extensions)]:
            basename = seq_fname.split(separator)[0] if separator else seq_fname[:seq_fname.rfind(".")]
            outpath = os.path.join(args.output_folder, args.project_folder, basename)
            os.makedirs(outpath, exist_ok=True)
            
            seq_path = os.path.join(input_path, seq_fname)
            # Parse sequence file
            genome = self.oIO.get_SeqDataset(seq_path, fname_separator = separator)
            contigs = genome.get_contigs()
            contig_titles = genome.get_contig_locations()
                
            if not genome or not contigs:
                raise ValueError(f"File {seq_fname} has wrong format or corrupted!")
                    
           
            # Parse GFF file    
            if master_gff is None:
                candidate_gff_files = [fn for fn in os.listdir(os.path.abspath(os.path.join(args.input_folder, args.project_folder))) 
                    if fn.upper().startswith(basename.upper()) and fn.upper().endswith(".GFF")]
                
                # Skip if pair sequence file has not been found
                if len(candidate_gff_files) == 0:
                    tools.alert(f"GFF file for sequence file {seq_fname} in folder {os.path.join(args.input_folder, args.project_folder)} has not been found!")
                    continue
                    
                   
                gff_fname = candidate_gff_files[0]
                gff_path = os.path.join(input_path, gff_fname)
                
            else:
                gff_fname = basename + ".gff"
                gff_path = os.path.join(outpath, gff_fname)
                gff_records = self.select_bin_gff(master_gff, genome.get_contig_titles())
                if len(gff_records) == 0:
                    tools.alert(f"GFF records for sequence file {seq_fname} have not been found!")
                    continue
                    
                self.oIO.save("\n".join(gff_records), gff_path)
            
            
            GFF_list = self.parse_gff_file(gff_path, contigs)
            
            
            # Groups of contigs
            group_collection = []
            # filtered contigs
            filtered_contigs = []
            try:
                filter_contigs = float(args.filter_chimeric_contigs)
                if filter_contigs > 0.05 or filter_contigs <= 0:
                    filter_contigs = False
            except:
                filter_contigs = True if args.filter_chimeric_contigs.upper()[0] == "Y" else False
            
            # Loop through motifs
            try:
                for motif in motifs:
                    options['-w'] = motif
                    tools.msg(f"Searching for motif {motif} in files {gff_fname} and {seq_fname}...")
                    # Run program
    
                    oMain = main.Interface(options, self.version, self.date_of_creation)
                    outputs = oMain.execute(genome, contigs, GFF_list)
                    
                    if not outputs:
                        continue
        
                    for output in outputs:
                        # groups = {'motif':motif, 'positive':[], 'negative':[], 'median':[], 'filtered_contigs':[], 'grouped':False}
                        groups = oMain.group_contigs(output=output, motif=motif, filter_contigs=filter_contigs)
                        filtered_contigs += groups['filtered_contigs']
            
                        if groups['grouped'] == False:
                            continue
                        
                        # Append grouped contigs
                        group_collection.append(groups)
                        
                        # Remove graphs from outputs if not needed
                        if args.save_graphs.upper() == "N" and 'svg' in output:
                            output['svg'] = None
                            
                        # Save output
                        oMain.save_output(output, path=outpath)
                        
                    # Clean TMP folder
                    oMain.clean_tmp_folder()
                    tools.msg("\tDone.")
            except:
                tools.msg("\tFailed!")
    
            grouped_contigs = self.group_contigs(contig_titles, group_collection, unresolved_contigs = args.unresolved_contigs)
            if len(grouped_contigs) == 0:
                tools.msg(f"Failed with {basename}!")
                continue
            
            seq = genome.get_sequence()
            for i in range(len(grouped_contigs)):
                contig_group = grouped_contigs[i]
                outfile_name = f"{basename}.{i + 1}.fasta" if len(grouped_contigs) > 1 else f"{basename}.fasta"
                group = []
                for j in range(len(contig_group)):
                    contig_location = contig_group[j]
                    contig = genome.get_contig_by_location(contig_location)
                    if contig is None:
                        sys.exit(f"ERROR: contig [{contig_location}] was not found!")
                    contig_title = str(contig.qualifiers['name'][0] if 'name' in contig.qualifiers else
                        (contig.qualifiers['id'][0] if 'id' in contig.qualifiers else f"contig.{i + 1}.{j + 1}"))
                        
                    # Filter contigs
                    if contig_title in filtered_contigs:
                        continue
                    group.append(f">{contig_title}\n{seq[contig.location.start:contig.location.end]}")
                # Save non-empty group to FASTA file
                if len(group):
                    self.oIO.save(data = "\n".join(group),
                        path = os.path.join(outpath, outfile_name),
                        data_format = "text"
                    )
                    
            # Save filtered contigs
            if filtered_contigs:
                self.oIO.save(data = "\n".join(filtered_contigs),
                    path = os.path.join(outpath, "filtered_contigs.txt"),
                    data_format = "text"
                )
 
