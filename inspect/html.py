import os, glob

from bamboo.html import png2base64str
from pandda.html import PANDDA_HTML_ENV
from pandda.constants import PanddaHtmlFilenames

def write_inspect_html(top_dir, inspector):

    # Get template to be filled in
    template = PANDDA_HTML_ENV.get_template('pandda_summary.html')
    # Output directory (for relative symlinks)
    out_dir  = os.path.join(top_dir, 'analyses', 'html_summaries')
    # Output file
    out_file = os.path.join(out_dir, PanddaHtmlFilenames.inspect_html)

    all_data = inspector.log_table
    len_data = len(all_data.index)

    #####################################################
    #                 Organise Data                     #
    #####################################################

    # Pandda Analysis outputs
    num_blobs = len_data
    num_sites = len(set(all_data['site_idx']))
    num_datasets = len(set(zip(*all_data.index)[0]))

    # Datasets Inspected/Modelled/Empty
    num_fitted = sum(all_data['Ligand Placed'])
    num_viewed = sum(all_data['Viewed'])
    num_empty  = num_viewed - num_fitted
    num_unviewed = len_data - num_viewed

    # Interesting unfitted (bookmarked)
    num_interesting_unfitted = sum(all_data["Interesting"][all_data['Ligand Placed']==False])

    # Number of Confident Ligand Models
    hgh_conf = sum(all_data["Ligand Placed"][all_data['Ligand Confidence']=='High'])
    med_conf = sum(all_data["Ligand Placed"][all_data['Ligand Confidence']=='Medium'])
    low_conf = sum(all_data["Ligand Placed"][all_data['Ligand Confidence']=='Low'])

    # Number of datasets with hits
    try:    num_d_hit = len(set(zip(*all_data.index[all_data['Ligand Placed'] == True])[0]))
    except: num_d_hit = 0
    # Number of sites with hits
    try:    num_s_hit = len(set(all_data[all_data['Ligand Placed'] == True]['site_idx']))
    except: num_s_hit = 0

    # Site Hits
    site_hits = []
    for site_idx in range(1, num_sites+1):
        site_placed = sum(all_data["Ligand Placed"][all_data['site_idx']==site_idx])
        if site_placed > 0: site_hits.append((site_idx, site_placed))

    #####################################################
    #               Prepare html data                   #
    #####################################################

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'PANDDA Inspect Summary'
    output_data['title'] = 'PANDDA Inspect Summary'
    output_data['introduction'] = 'Summary of Inspection of Datasets'
    # ===========================================================>
    # Header Images
    output_data['top_images'] = []
    if os.path.exists(os.path.join(out_dir, PanddaHtmlFilenames.pymol_sites_png_1)):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=os.path.join(out_dir, PanddaHtmlFilenames.pymol_sites_png_1))),
                                           'title': 'Identified Sites (Front)' })
    if os.path.exists(os.path.join(out_dir, PanddaHtmlFilenames.pymol_sites_png_2)):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=os.path.join(out_dir, PanddaHtmlFilenames.pymol_sites_png_2))),
                                           'title': 'Identified Sites (Back)' })
    for i_png, png in enumerate(sorted(glob.glob(os.path.join(out_dir, PanddaHtmlFilenames.inspect_site_graph_mult).format('*')))):
        output_data['top_images'].append({ 'path': 'data:image/png;base64,{}'.format(png2base64str(path=png)),
                                           'title': 'Identified Site Summary ({})'.format(i_png+1) })
    # ===========================================================>
    # Summary Bar
    output_data['summary_bar'] = []
    output_data['summary_bar'].append({'width':'4', 'colour':'success', 'text':'Datasets w. ligands: {} (of {})'.format(num_d_hit, num_datasets)})
    output_data['summary_bar'].append({'width':'4', 'colour':'success', 'text':'Sites w. ligands: {} (of {})'.format(num_s_hit, num_sites)})
    output_data['summary_bar'].append({'width':'4', 'colour':'warning', 'text':'Marked as interesting (unfitted): {}'.format(num_interesting_unfitted)})

    output_data['summary_bar'].append({'width':'3', 'colour':'info',    'text':'Total number of events: {}'.format(num_blobs)})
    output_data['summary_bar'].append({'width':'3', 'colour':'success', 'text':'High confidence hits:   {}'.format(hgh_conf)})
    output_data['summary_bar'].append({'width':'3', 'colour':'warning', 'text':'Medium confidence hits: {}'.format(med_conf)})
    output_data['summary_bar'].append({'width':'3', 'colour':'danger',  'text':'Low confidence hits:    {}'.format(low_conf)})
    # ===========================================================>
    # Progress Bars
    output_data['progress_bar'] = []
    # Fitting Progress
    output_data['progress_bar'].append({'title':'Fitting Progress', 'data':[]})
    output_data['progress_bar'][-1]['data'].append({'text':'Fitted - {} Events'.format(num_fitted),          'colour':'success', 'size':100.0*num_fitted/len_data})
    output_data['progress_bar'][-1]['data'].append({'text':'Unviewed - {} Events'.format(num_unviewed),      'colour':'warning', 'size':100.0*num_unviewed/len_data})
    output_data['progress_bar'][-1]['data'].append({'text':'No Ligand Fitted - {} Events'.format(num_empty), 'colour':'danger',  'size':100.0*num_empty/len_data})
    # Ligand Distribution
    if num_fitted > 0:
        import math
        weighting = sum([math.log1p(v[1]) for v in site_hits])
        output_data['progress_bar'].append({'title':'Identified Ligands by Site', 'data':[]})
        for i_block, (site_idx, n_hits) in enumerate(sorted(site_hits, key=lambda x: x[1], reverse=True)):
            colour = ('info','default')[i_block%2]
            output_data['progress_bar'][-1]['data'].append({'text':('S{}: {} hit'+'s'*bool(n_hits-1)).format(site_idx, n_hits), 'colour':colour, 'size':100.0*math.log1p(n_hits)/weighting})
    # ===========================================================>
    # Panels
    output_data['top_panels'] = []
    # Site summaries
    for site_idx, site_info in inspector.site_table.iterrows():
        if (site_info['Name'], site_info['Comment']) != ('None','None'):
            output_data['top_panels'].append({  'width':'6','color':'primary',
                                                'title':'<strong>Site {}</strong>: {}'.format(site_idx, site_info['Name']),
                                                'text':'<p><strong>Comment:</strong> '+site_info['Comment']+'</p>'  })
    if num_interesting_unfitted:
        output_data['top_panels'].append({  'width':'12','color':'primary',
                                            'title':'<strong>Interesting Datasets with No Ligand Placed</strong>',
                                            'text':'<p><strong>Datasets:</strong> '+', '.join(['<strong>{}</strong> (Event {})'.format(*i) for i in all_data[(all_data['Ligand Placed']==False)&(all_data['Interesting']==True)].index])+'</p>'  })
    # ===========================================================>
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = ['Dataset','Viewed','Interesting','Lig. Placed','Event','Site','1 - BDC','Z-Peak','Map Res.','Map Unc.','Confidence','Comment','']
    output_data['table']['rows'] = []
    # Add the datasets as rows
    for i_d in range(len(all_data.index)):

        # colour choices - 'success', 'muted', 'danger'
        # icon choices   - 'ok', 'flag', 'remove'

        d_data = all_data.iloc[i_d].to_dict()
        d_tag, d_event = all_data.index[i_d]

        columns = []

        if d_data['Viewed']:        columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Viewed']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Viewed']})

        if d_data['Interesting']:   columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Interesting']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Interesting']})

        if d_data['Ligand Placed']: columns.append({'colour':'success', 'icon':'ok',     'message':d_data['Ligand Placed']})
        else:                       columns.append({'colour':'danger',  'icon':'remove', 'message':d_data['Ligand Placed']})

        columns.append({'message':d_event})
        columns.append({'message':d_data['site_idx']})
        columns.append({'message':round(d_data['1-BDC'],3)})
        columns.append({'message':round(d_data['z_peak'],3)})
        columns.append({'message':d_data['analysed_resolution']})
        columns.append({'message':round(d_data['map_uncertainty'],3)})

        columns.append({'message':d_data['Ligand Confidence']})
        columns.append({'message':d_data['Comment']})

        row_message = 'Hit' if d_data['Ligand Placed'] else \
                      ''
        row_colour  = 'success' if d_data['Ligand Placed'] else \
                      'danger' if not d_data['Viewed'] else \
                      'info'

        output_data['table']['rows'].append({'heading' : d_tag,
                                             'colour'  : row_colour,
                                             'message' : row_message,
                                             'columns' : columns})

    with open(out_file, 'w') as out_html:
        out_html.write(template.render(output_data))
