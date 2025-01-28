#!/usr/bin/env python

from multiqc.modules.base_module import BaseMultiqcModule
import logging
import re
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger('multiqc.modules.agat')

class MultiqcModule(BaseMultiqcModule):
    print("AGAT module directory:", os.path.dirname(__file__))
    def __init__(self):
        super().__init__(
            name='agat',  # Change to lowercase
            anchor='agat',
            href="https://github.com/NBISweden/AGAT",
            info="is a suite of tools to handle gene annotations in any GTF/GFF format.",
            target="agat",
        )
        
        # Add debug prints
        print("AGAT module initialized")
        
        self.agat_data = dict()
        
        for f in self.find_log_files('agat'):
            parsed_data = self.parse_agat_stats(f['f'])
            if parsed_data:
                self.agat_data[f['s_name']] = parsed_data
                self.add_data_source(f)
            
        self.agat_data = self.ignore_samples(self.agat_data)
        
        if len(self.agat_data) == 0:
            raise UserWarning
            
        log.info(f"Found {len(self.agat_data)} reports")
        
        self.write_data_file(self.agat_data, 'multiqc_agat')
        self.add_general_stats()
        self.add_feature_type_plots()
    
    def parse_agat_stats(self, f):
        feature_types = ['region', 'mrna', 'ncrna', 'rrna', 'tmrna', 'trna']
        parsed_data = {feat_type: {} for feat_type in feature_types}
        current_section = None
        
        # Convert to TSV format
        lines = []
        for line in f.split('\n'):
            if line.startswith('----'):
                section_match = re.search(r'-+ ([a-z]+) -+', line)
                if section_match:
                    current_section = section_match.group(1)
                continue
            
            if current_section and line.strip():
                # Convert multiple spaces to tab
                tabbed_line = re.sub(r'\s{2,}', '\t', line.strip())
                if current_section in feature_types:
                    lines.append(f"{current_section}\t{tabbed_line}")
        
        return '\n'.join(lines)
    
    def add_general_stats(self):
        headers = {}
        metrics = [
            'Number of gene',
            'Number gene overlapping',
            'Total gene length (bp)',
            'median gene length (bp)',
            'Longest gene (bp)',
            'Shortest gene (bp)'
        ]
        
        for feat_type in ['mrna', 'ncrna', 'rrna', 'tmrna', 'trna']:
            for metric in metrics:
                header_key = f'{feat_type}_{metric}'
                display_name = metric.replace('(bp)', '').replace('Number of', '').strip()
                headers[header_key] = {
                    'title': f'{feat_type.upper()} {display_name}',
                    'description': f'{display_name} for {feat_type.upper()} genes',
                    'format': '{:,.0f}',
                    'suffix': ' bp' if 'bp' in metric else None,
                    'scale': 'RdYlBu' if 'overlapping' in metric else 'Blues'
                }
        
        data = {}
        for sample in self.agat_data:
            data[sample] = {}
            for feat_type in ['mrna', 'ncrna', 'rrna', 'tmrna', 'trna']:
                for metric in metrics:
                    data[sample][f'{feat_type}_{metric}'] = self.agat_data[sample][feat_type].get(metric, 0)
            
        self.general_stats_addcols(data, headers)
    
    def add_feature_type_plots(self):
        # Gene counts and overlaps plot
        for metric in ['Number of gene', 'Number gene overlapping']:
            counts_data = {}
            for sample in self.agat_data:
                counts_data[sample] = {}
                for feat_type in ['mrna', 'ncrna', 'rrna', 'tmrna', 'trna']:
                    counts_data[sample][feat_type.upper()] = self.agat_data[sample][feat_type].get(metric, 0)
            
            plot_config = {
                'id': f'agat_{metric.lower().replace(" ", "_")}',
                'title': f'AGAT: {metric} by Type',
                'ylab': 'Count',
            }
            
            self.add_section(
                name=f'{metric} by Type',
                anchor=f'agat_{metric.lower().replace(" ", "_")}',
                plot=bargraph.plot(counts_data, plot_config)
            )
        
        # Length statistics plots
        length_metrics = [
            'Total gene length (bp)',
            'median gene length (bp)',
            'Longest gene (bp)',
            'Shortest gene (bp)'
        ]
        
        for metric in length_metrics:
            length_data = {}
            for sample in self.agat_data:
                length_data[sample] = {}
                for feat_type in ['mrna', 'ncrna', 'rrna', 'tmrna', 'trna']:
                    length_data[sample][feat_type.upper()] = self.agat_data[sample][feat_type].get(metric, 0)
            
            metric_name = metric.replace('gene length (bp)', '').replace('gene (bp)', '').strip()
            plot_config = {
                'id': f'agat_gene_length_{metric_name.lower().replace(" ", "_")}',
                'title': f'AGAT: {metric_name} Gene Length by Type',
                'ylab': 'Length (bp)',
            }
            
            self.add_section(
                name=f'{metric_name} Gene Length',
                anchor=f'agat_gene_length_{metric_name.lower().replace(" ", "_")}',
                plot=bargraph.plot(length_data, plot_config)
            )
        
        # Region statistics
        if any('region' in d and d['region'] for d in self.agat_data.values()):
            region_data = {}
            region_metrics = [
                'Number of region',
                'Total region length (bp)',
                'mean region length (bp)',
                'median region length (bp)',
                'Longest region (bp)',
                'Shortest region (bp)'
            ]
            
            for sample in self.agat_data:
                region_data[sample] = {}
                for metric in region_metrics:
                    region_data[sample][metric] = self.agat_data[sample]['region'].get(metric, 0)
            
            region_config = {
                'id': 'agat_region_stats',
                'table_title': 'AGAT: Region Statistics',
                'sortRows': False,
                'col1_header': 'Sample'
            }
            
            self.add_section(
                name='Region Statistics',
                anchor='agat_region_stats',
                plot=table.plot(region_data, region_config)
            )