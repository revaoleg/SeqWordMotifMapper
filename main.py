import sys, os, re, time, copy
sys.path.append("lib")
import seq_io, tools, motifs, atlas, GI_finder, cui
from svg import SVG_dotplot

###############################################################################
# Command line interface
class Interface:
    def __init__(self, options=None, version="0", date_of_creation="<unknown>"):
        self.oValidator = cui.Validator()
        self.oIO = seq_io.IO()
        self.cwd = "."
        self.version = version
        self.date_of_creation = date_of_creation
            
        self.filtered_loci = []
        self.strand = 0
        self.oCompletedTasks = None
        
        if __name__ == "__main__":
            self.cwd = ".."

        self.options = cui.get_options()
        
        if options:
            self.options.update(options)
            
            if not self.options['-x']:
                self.options['-x'] = os.path.join(self.cwd,"lib","bin")
            if not self.options['-tmp']:
                self.options['-tmp'] = os.path.join(self.options['-x'],"tmp")
                
        self.input_folder = os.path.join(self.options['-u'], self.options['-d'])
        self.output_folder = os.path.join(self.options['-o'], self.options['-d'])
            
    def show_menu(self):
        oMenu = cui.Menu(self.options, version = self.version, date_of_creation = self.date_of_creation)
        self.options = oMenu.show()
        valid = self.oValidator.validate(self.options, echo=True)
        if valid:
            self.execute()
            
    def _parse_loci(self,fname):
        if not os.path.exists(fname):
            return []
        try:
            loci = self.oIO.read(fname,"text",inlist=True)
        except:
            return []

        if not loci:
            return []
        loci = [item.strip() for item in loci if item.strip() != ""]
        if not loci:
            return []
        success = False
        for symbol in ("-",".."," "):
            if loci[0].find(symbol) > 0:
                success = True
                break
        if not success:
            return loci
        try:
            loci = [[int(v) for v in item.split(symbol)] for item in loci if len(item.split(symbol)) == 2]
        except:
            return []
        return loci
    
    def save_options(self):
        self.oIO.save("\n".join(map(lambda item: "%s\t%s" % (item[0],str(item[1])), self.options.items())),os.path.join("lib","info"),"text")
        
    # Execute selected program
    def execute(self, genome, contigs, GFF_list):
        valid = self.oValidator.validate(self.options, echo=False)
        if not valid:
            tools.msg("Impropper option setting!")
            return None, None

        # Set input folder
        self.input_folder = os.path.join(self.options['-u'],self.options['-d'])

        # Parse regions to filter
        filtered_loci = []
        if self.options["-m"] and os.path.exists(os.path.join(self.input_folder,self.options['-m'])):
            try:
                filtered_loci = list(map(lambda ls: [int(ls[0]),int(ls[1])], 
                    self.oIO.read(os.path.join(self.input_folder,self.options['-m']),"text",inlist=True,separator="-",strip_symbol=" ")))
            except:
                raise ValueError(f"File {self.options['-m']} with filtered region boundaries does not exists or corrupted!")
        
        #### EXECUTION
        # Set strand
        self.strand = 1 if str(self.options['-std']).upper() in ("LEADING","DIR","1") else \
            (-1 if str(self.options['-std']).upper() in ("LAGGING","REV","-1") else 0)

        # Circular map
        if self.options['-mm'][0].upper() == "Y":
            return self.create_circular_motif_map(GFF_list=copy.deepcopy(GFF_list), genome_list=[genome], filtered_loci_list=[filtered_loci])

        # Stand-alone statistics panel
        if self.options['-mm'][0].upper() == "N" and self.options['-dp'][0].upper() == "N" and self.options['-sp'][0].upper() == "Y":
            return self.create_statplot(GFF_list=copy.deepcopy(GFF_list), genome_list=[genome], filtered_loci_list=[filtered_loci])
        
    def save_outputs(self, output_list, path = None):
        # Save output files
        for output in output_list:
            self.save_output(output, path)
            
    def save_output(self, output, path = None):
        # Set output folder
        if path is not None:
            output_folder = path
        else:
            output_folder = self.output_folder
        if output_folder is None:
            tools.msg("ERROR: output folder is not specified!")
            return   
        
        # Create output folder if does not exist
        os.makedirs(output_folder, exist_ok=True)
        
        # Save graphic file
        if output['svg']:
            if not output['svg output file']:
                output['svg output file'] = tools.ask_filename('SVG output file')
            if output['svg output file']:
                grap_output_path = os.path.join(output_folder,output['svg output file'])
                output['svg'].save(output_file=grap_output_path, output_format=self.options['-ogf'], dpi=int(self.options['-dpi']))
        # Save text output
        if output['report']:
            if not output['report output file']:
                output['report output file'] = tools.ask_filename('report TXT file')
            if output['report output file']:
                report_output_path = os.path.join(output_folder,output['report output file'])
                self.oIO.save(output['report'],report_output_path,"text")
                    
    def group_contigs(self, output, motif="", p_cutoff=0.05, filter_contigs=False):
        if 'report' not in output:
            sys.exit("Wrong format of output. Key 'report' is expected!")
            
        def is_distribution_Biased(data, p_cutoff):
            start = data.find("p-value") + 7
            stop = data.find(";")
            try:
                p = abs(float(data[start:stop]))
            except:
                return False
            
            if p <= p_cutoff:
                return True
            return False
            
        def parse_line(line):
            parts = line.strip().split(" ")
            try:
                title = parts[0][1:-1]  # start..end
                z = float(parts[2])
                p = float(parts[-1])
                return title, z, p
            except:
                sys.exit(f"Wrong line format: {parts}!")
                
        def is_distribution_biased(title, modified_sites, p_cutoff=None):
            coords = title.split("..") if title.find("..") > -1 else title.split("-")
            start, stop = [int(v) for v in coords]
            sites = [location for location in modified_sites if location >= start and location <= stop]
            if p_cutoff is not None:
                return tools.chi2_biased_distribution(sites, p_cutoff)
            return tools.median_biased_distribution(sites)
        
        # Group holder    
        groups = {'motif':motif, 'positive':[], 'negative':[], 'median':[], "filtered_contigs":[], 'grouped':False, 'distribution':''}
        data = output['report']
        modified_sites = []
        if filter_contigs:
            modified_sites = data[data.find("Motif\tStrand\tLocation"):].split("\n")[1:]
            # Store locations of modified sites
            # Motif	Strand	Location	Site	Match	Description	Score	Coverage	Site location relative to TSC	Annotation
            # CTAG	REV	132210..132213	132212	Yes	modified_base	37	4	n/a
            modified_sites = [int(line.split("\t")[3]) for line in modified_sites if len(line.split("\t")) == 9]
        
        # Locate contig distribution data
        block_start = data.find("Distribution across contigs:")
        distribution = "biased"
        block_end = data.find(f"The hypothesis of {distribution} distribution")
        distribution = "uniform"
        if block_end == -1:
            block_end = data.find("The hypothesis of {distribution} distribution")
        if block_start == -1 or block_end == -1:
            return groups
            
        groups['grouped'] = is_distribution_Biased(data[block_end:], p_cutoff)
        if groups['grouped'] == False:
            return groups
        # Assign groups
        groups['distribution'] = distribution
        block = data[block_start:block_end].split("\n")
        
        for line in block:
            if line.strip().startswith("["):
                title, z, p = parse_line(line)
                if filter_contigs:
                    if isinstance(filter_contigs, float):
                        p_cutoff = filter_contigs
                    else:
                        p_cutoff = None
                    if is_distribution_biased(title, modified_sites, p_cutoff):
                        groups['filtered_contigs'].append(title)
                if p > p_cutoff:
                    groups['median'].append(title)
                elif z < 0:
                    groups['negative'].append(title)
                else:
                    groups['positive'].append(title)
        return groups
    
    # Prepare string description of motifs used for modified base filtering
    def get_motif_description(self, graph_name):
        if graph_name == "cmap":
            return self.options['-w']
        else:
            nucleotides = [s.strip().upper() for s in self.options['-dpn'].split(',')]
            mtypes = [s.strip() for s in self.options['-dpm'].split(',')]
            return (",".join(nucleotides) if nucleotides else ",".join(mtypes)) + f" {self.options['-dpf']}"
            
    # Prepare outfile name
    def get_outfilename(self, gff_filename, graph_name):
        # Prepare prefix
        def _prefix():
            prx = ""
            if self.options['-s'][0].upper() == "M":
                prx = "m"
            if self.options['-f'][0].upper() == "U":
                prx = "u"+prx
            return prx
        generic_filename =  gff_filename[:-gff_filename.rfind(".")] if not self.options['-ft'] else self.options['-ft']
        generic_filename = tools.check_file_name_symbols(generic_filename)
        motif_description = self.get_motif_description(graph_name)
        return f"{_prefix() if graph_name=='cmap' else ''}{generic_filename}.{graph_name}_{motif_description}"

    # Filter GFF records according to program run settings
    def filter_GFF_records(self, gff_records, cutoff=0, filtered_regions=[], nucleotides=[], mtypes=[], motifs=[], flg_exclude_strand='OFF'):
        # Filter GFF records by score cutoff values
        if cutoff:
            gff_records = [record for record in gff_records if float(record['score']) >= cutoff]
            
        # Filter sites from excluded regions
        if not filtered_regions:
            filtered_regions = self.filtered_loci
        if filtered_regions:
            gff_records = tools.filter_regions(sites=gff_records, regions=filtered_regions)
                
        # Filter GFF records by modified nucleotides
        if nucleotides and nucleotides[0]:
            gff_records = list(filter(lambda dot: dot['nucleotide'] in nucleotides, gff_records))
            
        # Filter GFF records by types of methylation
        if mtypes and mtypes[0]:
            gff_records = list(filter(lambda dot: dot['modification'] in mtypes, gff_records))
        
        # Exclude methylation on the reverse-complement strand
        if flg_exclude_strand.upper() in ('LEADING', '+', '1'):
            gff_records = [d for d in gff_records if d['strand'] == "+"]
        if flg_exclude_strand.upper() in ('LAGGING', '-', '-1'):
            gff_records = [d for d in gff_records if d['strand'] == "-"]
        
        # Filter GFF records by motifs
        '''
        # Filtering is needed only for dotplot and standalong statplot, as circular plot does it internaly
        if motifs and (self.options['-dp'][0].upper() == "Y" or 
            (self.options['-dp'][0].upper() == "N" and self.options['-mm'][0].upper() == "N" and self.options['-sp'][0].upper() == "Y")):
                gff_records = tools.filter_motifs(records = gff_records, motifs = [motif.strip() for motif in motifs.split(";")])
        '''
        if motifs:
            gff_records = tools.filter_motifs(records = gff_records, motifs = [motif.strip() for motif in motifs.split(";")])
        return gff_records
    
    # Create list of motif object for further verification
    def set_MotifObjLs(self, motif_string,    # Motifs = 'GATC,2,-2; AGNCT,1,-1'
        reference, binpath, tmppath, 
        score_cutoff=20, promoter_length = 75, context_mismatches = 0, motif_mismatch = 0, 
        motifs_or_sites='sites', modified_or_unmodified_sites="M"):
        MotifObjLs = []
        # Circle through motifs
        for motif in [mf.strip() for mf in motif_string.split(";")]:
            # Ignore negative motifs, they will be filtered out later 
            if motif.startswith("-"):
                continue
            # Parse motif setting
            word, direct_locations, reverse_locations = self.parse_motif(motif)
            # Check option for searching either methlated or unmodified sites
            find_modified_sites = False
            if self.options['-f'] == "M":
                find_modified_sites = True
            # Create a motif object
            oMotifObj = motifs.Motif(reference=reference, binpath=binpath, tmppath=tmppath, motif=motif,
                word=word, motifs_or_sites=motifs_or_sites, modbase_location=direct_locations, reverse_modbase_location=reverse_locations,
                score_cutoff=score_cutoff, promoter_length=promoter_length, context_mismatches=context_mismatches, strand=self.strand,
                motif_mismatch=motif_mismatch, filtered_loci=self.filtered_loci, modified_or_unmodified_sites=modified_or_unmodified_sites,
                max_sites_for_verification=self.options['-r'])
            # Add motif object to the list
            MotifObjLs.append(oMotifObj)
        return MotifObjLs
            
    def create_circular_motif_map(self, GFF_list, genome_list, filtered_loci_list=[]):
        # Cycle across genomes
        outputs = []
        for g in range(len(genome_list)):
            genome = genome_list[g]
            self.filtered_loci = filtered_loci_list[g]

            # Create a list of motifs
            if os.path.exists(os.path.join(self.options['-u'], self.options['-d'], self.options['-w'])):
                oMotifObjLs = _parse_motif_list()
            else:
                oMotifObjLs = self.set_MotifObjLs(motif_string=self.options['-w'],    # Motifs = 'GATC,2,-2; AGNCT,1,-1'
                    reference=genome.get_description(), binpath=self.options['-x'], tmppath=self.options['-tmp'], 
                    score_cutoff=self.options['-c'], promoter_length=self.options['-p'], context_mismatches=self.options['-n'], motif_mismatch=self.options['-z'], 
                    motifs_or_sites=self.options['-s'], modified_or_unmodified_sites=self.options['-f'].upper()[0])
    
            for oMotifObj in oMotifObjLs:
                if oMotifObj == None:
                    continue
                for gff in GFF_list:
                    
                    if isinstance(gff, dict):
                        # Filter modified sites according to user specuified parameters
                        gff['Body'] = self.filter_GFF_records(gff_records=gff['Body'], cutoff=float(self.options['-c']), motifs=self.options['-w'],
                            flg_exclude_strand=self.options['-std'])
                        # Check if the number of filtered sites lower than the maximum number of sites for verification
                        if int(self.options['-r']) and len(gff['Body']) > int(self.options['-r']):
                            raise ValueError(f"Number of filtered sites {len(gff['Body'])} is bigger than the maximum number\n" +
                            f"of sites for verification - {self.options['-r']}!\n" +
                            f"To fix the problem, option '-r' can be increased or set to zero. However, it may cause problems with program run.\n" +
                            f"A better approach is to use more stringent filter settings or more complex motifs to reduce the number of selected sites.")
                    
                    oMotifObj.reset()
                    oMotifObj.execute(gff, genome.get_genome())

                    outfile = self.get_outfilename(gff_filename=gff['input file'], graph_name="cmap")
                    if self.options['-s'] == 'motifs':
                        msg_found = "%d/%d" % (oMotifObj.num_found_motifs,oMotifObj.num_expected_sites)
                    elif self.options['-f'].upper()[0] == "M":
                        msg_found = "%d/%d" % (oMotifObj.num_found_sites,oMotifObj.num_expected_sites)
                    elif self.options['-f'].upper()[0] == "U":
                        msg_found = "%d/%d" % (oMotifObj.num_found_sites,oMotifObj.num_expected_sites)
                    else:
                        msg_found = ""
                    
                    oAtlas = atlas.Main(seqfile = genome,
                        title = oMotifObj.get_word(),
                        modbases = oMotifObj.get_modbase_dict(),
                        loci = self.filtered_loci,
                        window_length = self.options['-wl'],
                        window_step = self.options['-ws'],
                        graph_format = self.options['-ogf'].upper(),
                        graph_title = self.options['-cmt'],
                        tmp_folder = self.options['-tmp'])
                        
                    oSVG = oAtlas.svg([oMotifObj.get_modbase_location(),oMotifObj.get_rev_modbase_location(),msg_found])
                    
                    # Add statplot if requested
                    statplot_report = ""
                    if self.options['-sp'].upper()[0] == "Y":
                        sites = oMotifObj.get_modbase_list()
                        
                        # Check if the number of sites is sufficient for statistics
                        if len(sites) >= 10:
                            # Estimate top indend for statplot
                            top_indend = 0
                            if oSVG:
                                top_indend = oSVG.get_height()
                            # Receive a list of modified nucleotides/motifs in the format [[start,end,strand,{data}],...]
                            # Create statplot SVG
                            oStatplotSVG, statplot_report, oTasks = self.append_statplot(genome=genome, sites=sites, top_indend=top_indend, 
                                motif_setting=self.options['-w'], verified=True)
                            self.oCompletedTasks = oTasks
                            # Add statplot SVG to dotplot SVG
                            if oStatplotSVG:
                                if 'contigs' in oTasks and 'contig borders' in oTasks['contigs']:
                                    oSVG.format_contigs(oTasks['contigs']['contig borders'])    # [{title, start, end, z-score : (score, err, p, exp, found)}...
                                oSVG = oSVG.join(oStatplotSVG, self_tag="circular_map", self_transform="h_center", other_tag="stat_plot")
                                
                    outputs.append({'input file':gff['input file'], 'svg':oSVG, 'report':oMotifObj.tostring(text_to_insert=statplot_report),
                        'svg output file':outfile, 'report output file':outfile + ".txt"})
        return outputs
                        
    def create_dotplot(self, GFF_list, genome_list, filtered_loci_list=[]):
        # Processing GFF data
        def _GFF_execute(gff, genome, width, height, title="", motif_setting="", verified=False):
            dots = gff['Body']
            report = ""
            # Set maximal X (coverage) and Y (score) values
            try:
                if not width:
                    width = max([int(dots[i]['data']['coverage']) for i in range(len(dots))])
                    width = int(str(width)[0]+("0"*(len(str(width))-1)))+10**(len(str(width))-1)
                if not height:
                    height = max([int(dots[i]['score']) for i in range(len(dots))])
                    height = int(str(height)[0]+("0"*(len(str(height))-1)))+10**(len(str(height))-1)
            except:
                tools.alert(f"File {gff['input file']} is either empty, or corrupted, or all modified nucleotides were filtered out!")
                return 
            
            # Generate DotPlot SVG file
            #svg_data = svg(data=dots, title=title, width=width, height=height)
            svg_data = SVG_dotplot(data=dots, title=title, graph_title=self.options['-dpt'], width=width, height=height)
            svg_data.set_svg()
            
            # Add statplot
            statplot_report = ""
            if self.options['-sp'].upper()[0] == "Y" and len(dots) > 10 and verified:
                # Estimate top indend for statplot
                top_indend = 0
                if svg_data:
                    top_indend = svg_data.get_height()
                # Receive a list of modified nucleotides/motifs in the format [[start,end,strand],...]
                sites = [[int(site['start']), int(site['end']), 1 if site['strand'] == '+' else -1, 
                    {'Annotation':site['data']['Annotation'], 'Location relative to TSS':site['data']['Location relative to TSS'], 'Site':site['site']}
                        if ('Annotation' in site['data'] and 'Location relative to TSS' in site['data'] and 'site' in site) else {}]
                    for site in dots]
                # Create statplot SVG
                oStatplotSVG, statplot_report, oTasks =  self.append_statplot(sites=sites, top_indend=top_indend, 
                    motif_setting=motif_setting, verified=verified)
                self.oCompletedTasks = oTasks
                # Add statplot SVG to dotplot SVG
                if oStatplotSVG:
                    svg_data += oStatplotSVG
            return svg_data, report + statplot_report
            
        # Select sought bases
        def _get_sougth_bases(nucleotides, mtypes, motifs):
            bases = [Nuc for Nuc in nucleotides]
            mt_bases = [mtype[1] for mtype in mtypes if len(mtype)> 1]
            if len(mt_bases):
                bases = mt_bases
            motif_items = [motif.strip() for motif in motifs if not motif.strip().startswith("-")]
            for item in motif_items:
                parts = [part.strip() for part in item.split(",")]
                if len(parts) > 1:
                    word = parts[0].upper()
                    points = [int(v) for v in parts[1:] if (v and int(v))]
                    bases += [word[p - 1] if p > 0 else word[p] for p in points]
            bases = tools.dereplicate(bases)
            return bases
        
        outputs = []
        for g in range(len(genome_list)):
            # Get genome object 
            genome = genome_list[g]
            # Set filtered loci
            self.filtered_loci = filtered_loci_list[g]
            for gff in GFF_list:
                # Set output file names
                INPUT_FILE = gff['input file']
                GENERIC_FILE_NAME = self.get_outfilename(gff_filename=INPUT_FILE, graph_name='dotplot')
                GRAPH_OUTPUT_FILE = GENERIC_FILE_NAME
                REPORT_FILE = GENERIC_FILE_NAME + ".txt"
                
                # Filter modified sites according to user specuified parameters
                nucleotides = [s.strip().upper() for s in self.options['-dpn'].split(',')]
                mtypes = [s.strip() for s in self.options['-dpm'].split(',')]
                motif_set = self.options['-dpf']
                gff['Body'] = self.filter_GFF_records(gff_records=gff['Body'], 
                    cutoff=float(self.options['-dpc']), nucleotides = nucleotides,
                    mtypes = mtypes, motifs = motif_set, flg_exclude_strand=self.options['-std'])

                # If statplot is requested, check whether the number of filtered sites lower than the maximum number of sites for verification
                flg_verify = False
                if self.options['-sp'][0].upper() == "Y":
                    flg_verify = True
                    if int(self.options['-r']) and len(gff['Body']) > int(self.options['-r']):
                        flg_verify = False
                        tools.alert(f"Number of identified sites, {len(gff['Body'])}, exeeds the verification limit {self.options['-r']}!")
                    if flg_verify:
                        oMotifObj = motifs.Motif(reference=genome['description'], binpath=self.options['-x'], tmppath=self.options['-tmp'],
                            promoter_length=self.options['-p'], context_mismatches=self.options['-n'], motif_mismatch=self.options['-z'],
                            filtered_loci=self.filtered_loci, strand=self.strand)
                        success = oMotifObj.execute(gff_data = os.path.join(self.input_folder,gff) if isinstance(gff, str) else gff,
                            gbk_data = os.path.join(self.input_folder,genome) if isinstance(genome, str) else genome)
                        if not success:
                            return
                        
                        gff['Body'] = oMotifObj.get_entries()
                        '''
                        # Check missmatches between expected and identified bases
                        sought_bases = _get_sougth_bases(nucleotides, mtypes, motif_set)
                        gff['Body'] = [site for site in gff['Body'] if site['nucleotide'].upper() in sought_bases]
                        '''
    
                # Check if any sites were selected after filtering
                if len(gff['Body']):
                    # Generate SVG code and report
                    oSVG, report = _GFF_execute(gff=gff, genome=genome,
                        title=self.options['-dpt'], width=self.options['-dpw'], height=self.options['-dps'], motif_setting=self.get_motif_description('dotplot'), verified=flg_verify
                    )
                    # Collect output data
                    outputs.append({'input file':gff['input file'], 'svg':oSVG, 'report':report,
                        'svg output file':GRAPH_OUTPUT_FILE, 'report output file':REPORT_FILE})
        return outputs
    
    # Create a stand-alone statistical plot
    def create_statplot(self, GFF_list, genome_list, filtered_loci_list=[]):
        outputs = []
        for g in range(len(genome_list)):
            # Get genome object 
            genome = genome_list[g]
            # Set filtered loci
            self.filtered_loci = filtered_loci_list[g]
            for gff in GFF_list:
                # Set output file names
                INPUT_FILE = gff['input file']
                GENERIC_FILE_NAME = self.get_outfilename(gff_filename=INPUT_FILE, graph_name='statplot')
                GRAPH_OUTPUT_FILE = GENERIC_FILE_NAME
                REPORT_FILE = GENERIC_FILE_NAME + ".txt"
                
                # Filter modified sites according to user specuified parameters
                nucleotides = [s.strip().upper() for s in self.options['-dpn'].split(',')]
                mtypes = [s.strip() for s in self.options['-dpm'].split(',')]
                # Filtering
                sites = self.filter_GFF_records(gff_records=gff['Body'], 
                    cutoff=float(self.options['-dpc']), nucleotides = nucleotides,
                    mtypes = mtypes, motifs = self.options['-dpf'], flg_exclude_strand=self.options['-std'])

                # Check if the number of filtered sites is sufficient for statistical analysis
                if len(sites) < 10:
                    raise ValueError(f"Number of modified sites is {len(sites)}. Fos statistical analysis, number of sites should be > 10.") 
                # Check if the number of filtered sites is lower than the maximum number of sites for verification
                flg_verify = True
                if int(self.options['-r']) and len(sites) > int(self.options['-r']):
                    flg_verify = False
                if flg_verify:
                    # Verification of sites
                    oMotifObj = motifs.Motif(reference=genome['description'], binpath=self.options['-x'], tmppath=self.options['-tmp'],
                        promoter_length=self.options['-p'], context_mismatches=self.options['-n'], motif_mismatch=self.options['-z'],
                        filtered_loci=self.filtered_loci, strand=self.strand)
                    success = oMotifObj.execute(gff_data = sites,
                        gbk_data = os.path.join(self.input_folder,genome) if isinstance(genome, str) else genome)
                    sites = oMotifObj.get_entries()

                # Receive a list of modified nucleotides/motifs in the format [[start,end,strand,{data}],...]
                sites = [[int(site['start']), int(site['end']), 1 if site['strand'] == '+' else -1, 
                    {'Annotation':site['data']['Annotation'], 'Location relative to TSS':site['data']['Location relative to TSS'], 'Site':site['site']}
                        if ('Annotation' in site['data'] and 'Location relative to TSS' in site['data'] and 'site' in site) else {}]
                    for site in sites]
                # Create statplot SVG
                
                PlotSVG, report, _ =  self.append_statplot(sites=sites, motif_setting=self.get_motif_description('statplot'), 
                    verified=flg_verify, graph_title=self.options['-spt'])

                # Collect output data
                outputs.append({'input file':gff['input file'], 'svg':PlotSVG, 'report':report,
                    'svg output file':GRAPH_OUTPUT_FILE, 'report output file':REPORT_FILE})
        return outputs
        
    # Append statistical plot
    def append_statplot(self, genome=None, sites=[], statplot_setting=None, top_indend=50, motif_setting="", verified=False, graph_title=""):
        report = ""
        if not statplot_setting:
            statplot_setting = self.generate_statplot_setting()
        #Generate DensityPlot SVG
        if statplot_setting:
            oDensityPlot = GI_finder.Interface(options=statplot_setting, completed_tasks=self.oCompletedTasks, motif_setting=motif_setting, 
                max_sites_for_verification=self.options['-r'], graph_title=graph_title)
            oDensityPlot.set_task("MOD", flg_set_as_immutable=True)
            # Set tasks
            if self.options['-tsk']:
                tasks = [s.strip().upper() for s in self.options['-tsk'].split(",")]
                for task in tasks:
                    oDensityPlot.set_task(task)
            # Pass input GBK file and locations of modified nucleotides
            GBK_files = statplot_setting['-g'] if genome is None else genome
            PlotSVG, report = oDensityPlot.execute(GBK_files = GBK_files, top_indend = top_indend, 
                modified_sites=[sites], verified=verified)
            return PlotSVG, report, oDensityPlot.getCompletedTasks()
        return None, report, None
        
    def generate_statplot_setting(self):
        oPlotSetting = GI_finder.PlotSettings()
        statplot_settings = oPlotSetting(default=True)
        if statplot_settings:
            # Set path to GBK file
            statplot_settings['-g'] = os.path.join(self.input_folder,self.options['-g'])
            # Leading/Lagging strands
            if str(self.options['-std']).upper() in ('LEADING', '1', 'DIR'):
                strand = "leading"
            elif str(self.options['-std']).upper() in ('LAGGING', '-1', 'REV'):
                strand = "lagging"
            else:
                strand = "off"
            statplot_settings['-std'] = strand
            # Promoter length
            statplot_settings['-p'] = self.options['-p']
            # Motif setting
            if self.options['-mm'].upper()[0] == "Y":
                statplot_settings['-w'] = self.options['-w']
        return statplot_settings
            
    #### GFF PROCESSING FUNCTIONS
    def parse_motif(self, motif):
        if not motif.strip():
            return "","",""
        values =  [s.strip() for s in motif.split(",")]
        try:
            word = values[0]
            values = [int(v) for v in values[1:]]
        except:
            raise ValueError(f"Wrong motif {motif}!")
        direct_locations = ",".join([str(v) for v in values if v > 0])
        reverse_locations = ",".join([str(len(word) + v + 1) for v in values if v < 0])
        return word, direct_locations if direct_locations else "", reverse_locations if reverse_locations else ""
    
    def parse_included_excluded_motifs(self, motifs):
        if not motifs:
            return [], []
        # Sort out the included and excluded motifs by sign '-' in front of the motifs, like -GATC,2,-2; AGCT,1,-1
        motifs_to_include = [motif.strip() for motif in motifs.split(";") if not motif.strip().startswith("-")]
        motifs_to_exclude = [motif[1:].strip() for motif in motifs.split(";") if motif.strip().startswith("-")]
        # Split each motif by commas to motifs and numbers of modified nucleotides
        motifs_to_include = [s.split(",") for s in motifs_to_include]
        motifs_to_exclude = [s.split(",") for s in motifs_to_exclude]
        # Convert modified nucleotide locations to integers
        motifs_to_include = [[ls[0].strip().upper()] + [int(v) for v in ls[1:]] for ls in motifs_to_include]
        motifs_to_exclude = [[ls[0].strip().upper()] + [int(v) for v in ls[1:]] for ls in motifs_to_exclude]    
        return motifs_to_include, motifs_to_exclude
        
    def clean_tmp_folder(self):
        tmp_folder = self.options["-tmp"]
        current_time = time.time()
    
        for fname in os.listdir(tmp_folder):
            file_path = os.path.join(tmp_folder, fname)
            # Check if the file is at least 24 hours old (86400 seconds)
            if os.path.isfile(file_path) and (current_time - os.path.getmtime(file_path)) > 86400:
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}: {e}")            
