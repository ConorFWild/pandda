import os

from scitbx.array_family import flex

from giant.html import GIANT_HTML_ENV
from bamboo.html import path2url

def write_html_summary(fname, atom_group_summaries):

    # Get template to be filled in
    template = GIANT_HTML_ENV.get_template('summary_page.html')

    # ===========================================================>
    # Construct the data object to populate the template
    output_data = {}
    output_data['header'] = 'Structure Summary'
    output_data['title'] = 'Structure Summary'

#    # ===========================================================>
#    # Progress Bars
#    output_data['progress_bar'] = []
#    output_data['progress_bar'].append({'title':'Dataset Summary', 'data':[]})
#    output_data['progress_bar'][0]['data'].append({'text':'Interesting',     'colour':'success', 'size':100.0*num_interesting/len_data})
#    output_data['progress_bar'][0]['data'].append({'text':'Not Interesting', 'colour':'info',    'size':100.0*num_not_interesting/len_data})
#    output_data['progress_bar'][0]['data'].append({'text':'Not Analysed',    'colour':'warning', 'size':100.0*num_not_analysed/len_data})
#    output_data['progress_bar'][0]['data'].append({'text':'Rejected',        'colour':'danger',  'size':100.0*num_rejected/len_data})

    # ===========================================================>
    # Tables
    output_data['table'] = {}
    output_data['table']['column_headings'] = [ 'Num Atoms',

                                                'Min BFactor',
                                                'Max BFactor',
                                                'Mean BFactor',
                                                'RMS BFactor',

                                                'Min BFactor-Z',
                                                'Max BFactor-Z',
                                                'Mean BFactor-Z',
                                                'RMS BFactor-Z',
                                            ]
    output_data['table']['rows'] = []

    for ag_sum in atom_group_summaries:
        # Extract the atom group info as a dictionary
        ed_sum = ag_sum.edstats_scores

        bfs      = ag_sum.b_factors
        bfs_stat = ag_sum.b_factors_statistics
        bfz      = ag_sum.b_factors_z
        bfz_stat = ag_sum.b_factors_z_statistics

        # colour choices - 'success', 'muted', 'danger'
        # icon choices   - 'ok', 'flag', 'remove'

        columns = []
        # Num Atoms
        columns.append({'colour':'default', 'message':round(bfs.size(),1)})
        # B-Factors
        columns.append({'colour':'default', 'message':round(bfs_stat.min,1)})
        columns.append({'colour':'default', 'message':round(bfs_stat.max,1)})
        columns.append({'colour':'default', 'message':round(bfs_stat.mean,1)})
        columns.append({'colour':'default', 'message':round((flex.sum_sq(bfs)/bfs.size())**0.5,1)})

        # B-Factors Z-Scores
        columns.append({'colour':'default', 'message':round(bfz_stat.min,2)})
        columns.append({'colour':'default', 'message':round(bfz_stat.max,2)})
        columns.append({'colour':'default', 'message':round(bfz_stat.mean,2)})
        if (flex.sum_sq(bfz)/bfz.size())**0.5 < 2.5:
            columns.append({'colour':'default', 'message':round((flex.sum_sq(bfz)/bfz.size())**0.5,1)})
        else:
            columns.append({'colour':'danger', 'icon':'remove', 'message':round((flex.sum_sq(bfz)/bfz.size())**0.5,1)})

        row_message = '' if (flex.sum_sq(bfz)/bfz.size())**0.5 < 2.5 else \
                      'Check'
        row_colour  = 'success' if row_message == '' else \
                      'danger'

        output_data['table']['rows'].append({'heading' : ag_sum.atom_group.id_str(),
                                             'message' : row_message,
                                             'colour'  : row_colour,
                                             'columns' : columns})

    with open(fname, 'w') as out_html:
        out_html.write(template.render(output_data))
