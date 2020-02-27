import pandas as pd

import scipy.cluster

from scitbx.array_family import flex

from bamboo.common import Info
from bamboo.stats.cluster import find_connected_groups, generate_group_idxs

from pandda_2.filter_clusters import Event, PointCluster


class CreateSitesTable:
    def __init__(self,
                 clustering_cutoff=1.732,
                 make_pymol_site_image_and_scripts=None,
                 ):
        self.clustering_cutoff = clustering_cutoff
        self.make_pymol_site_image_and_scripts = make_pymol_site_image_and_scripts

    def __call__(self,
                 events_table,
                 grid,
                 reference,
                 ):
        """Cluster events to sites and add information to the pandda tables"""

        events = events_table_to_events(events_table)

        if not events:
            print('No Events Found')
            return None

        site_list = cluster_events(events=events,
                                   cutoff=self.clustering_cutoff / grid.grid_spacing(),
                                   linkage='average',
                                   )

        site_list.sort(key=lambda s: (s.info.num_events,
                                      max([e.cluster.max
                                           for e
                                           in s.children
                                           ]
                                          )
                                      ),
                       reverse=True).renumber()

        # Add meta to the site list
        # TODO implement this function -- blank at the moment TODO
        [s.find_protein_context(hierarchy=reference.model.hierarchy)
         for s
         in site_list.children
         ]

        records = []

        for site in site_list.children:
            record = {"site_idx": site.id,
                      "centroid": site.info.centroid,
                      "approx_size": site.info.approx_size,
                      "num_events": len(site.children),
                      "nearest_residue_1": None,
                      "nearest_residue_2": None,
                      "nearest_residue_3": None,
                      "near_crystal_contacts": site.info.near_crystal_contacts,
                      "native_centroid": get_native_centroid(site,
                                                             reference,
                                                             grid,
                                                             )
                      }
            records.append(record)

        sites_table = pd.DataFrame(records)

        events_table["dtag"] = events_table["dtag"].astype(str)
        events_table["event_idx"] = events_table["event_idx"].astype(str)
        events_table_with_sites = events_table.set_index(["dtag", "event_idx"])

        print(events_table)

        for site in site_list.children:
            for event in site.children:
                print(event.id)
                events_table_with_sites.loc[(str(event.id[0]), str(event.id[1])), "site_idx"] = site.id

        return sites_table, events_table_with_sites

        #
        # # Update the pandda tables?
        # if update_tables:
        #     self.update_site_table(site_list=site_list, clear_table=True)
        #     self.update_event_table_site_info(events=events)
        #
        # # Generate output images and graphs?
        # if update_output:
        #     # Plot output graph of site list
        #     self.log('Deleting old images: ')
        #     delete_with_glob(glob_str=self.file_manager.get_file('analyse_site_graph_mult').format('*'))
        #     bar.multiple_bar_plot_over_several_images(
        #         f_template=self.file_manager.get_file('analyse_site_graph_mult'),
        #         plot_vals=[sorted([e.cluster.max for e in s.children], reverse=True) for s in site_list.children])
        #     # Create pictures of the sites on the protein
        #     if self.make_pymol_site_image_and_scripts:
        #         self.make_pymol_site_image_and_scripts(site_list=site_list, make_images=True)


def get_native_centroid(site,
                        reference,
                        grid,
                        ):
    centroid_cart = tuple(flex.double(site.info.centroid) * grid.grid_spacing())

    return tuple(reference.model.alignment.ref2nat(coordinates=[centroid_cart])[0])


def cluster_events(events,
                   cutoff=10,
                   linkage='average',
                   ):
    if len(events) == 1:
        return SiteList([Site(events=events, id=1).apply_parentage()])

    centroids = [e.cluster.centroid
                 for e
                 in events
                 ]
    cluster_ids = scipy.cluster.hierarchy.fclusterdata(X=centroids,
                                                       t=cutoff,
                                                       criterion='distance',
                                                       metric='euclidean',
                                                       method=linkage,
                                                       )
    cluster_ids = list(cluster_ids)

    sites = []
    for s_idx, e_idxs in generate_group_idxs(cluster_ids):
        assert s_idx > 0
        new_site = Site([events[i]
                         for i
                         in e_idxs
                         ],
                        id=s_idx).apply_parentage()
        sites.append(new_site)

    return SiteList(sites)


def update_site_table(self, site_list, clear_table=True):
    """Add site entries to the site table"""

    # Clear an existing table
    if clear_table:
        self.tables.site_info = pd.DataFrame(data=None,
                                             index=self.tables.site_info.index.reindex([])[0],
                                             columns=self.tables.site_info.columns)
    # Go through and update the site information
    for site in site_list.children:
        self.tables.site_info.loc[site.id, :] = None
        centroid_cart = tuple(flex.double(site.info.centroid) * self.grid.grid_spacing())
        self.tables.site_info.set_value(site.id, 'centroid', centroid_cart)
        self.tables.site_info.set_value(site.id, 'native_centroid', tuple(
            self.datasets.reference().model.alignment.ref2nat(coordinates=[centroid_cart])[0]))


# def update_event_table_site_info(self, events):
#     """Update the event table for pre-existing events"""
#     for e in events:
#         assert e.id, 'NO ID GIVEN: {!s}'.format(e.id)
#         assert e.parent, 'EVENT HAS NO PARENT: {!s}'.format(e.parent)
#         self.tables.event_info.set_value(e.id, 'site_idx', e.parent.id)

class Site(object):
    _attributes = ['centroid', 'num_events', 'approx_size',
                   'nearest_residues', 'near_crystal_contacts']

    def __init__(self, events=None, id=None, info=None):
        """Class to hold information about an identified site (collection of events)"""
        # List of Events
        self.children = []
        self.parent = None
        self.id = id
        # Add Meta to the object
        if info:
            assert isinstance(info, Info)
            for a in self._attributes: assert hasattr(info, a)
            self.info = info
        else:
            self.info = Info(self._attributes)

        # Add Events
        if events: self.add_events(events=events, update=True)

    def add_events(self, events, update=True):
        if isinstance(events, list):
            assert isinstance(events[0], Event), 'Added events must be of class Event'
            self.children.extend(events)
        else:
            assert isinstance(events, Event), 'Added events must be of class Event'
            self.children.append(events)
        if update: self.update()
        return self

    def apply_parentage(self):
        """Set the site as the parents of the events"""
        for e in self.children: e.parent = self
        return self

    def update(self):
        centroids = flex.vec3_double([e.cluster.centroid
                                      for e
                                      in self.children
                                      ]
                                     )
        self.info.centroid = centroids.mean()
        self.info.approx_size = (centroids.min(), centroids.max())
        self.info.num_events = len(self.children)
        return self

    def summary(self):
        out = []
        out.append('Site Summary')
        out.append('{!s} Events'.format(self.info.num_events))
        out.append('Centroid: {!s}'.format(self.info.centroid))
        out.append('Approx Range: {!s}'.format(self.info.approx_size))
        for s in self.children:
            out.append('-----------------------')
            out.append(s.summary())
        return '\n'.join(out)

    def sort(self, key, reverse=True):
        self.children.sort(key=key, reverse=reverse)
        return self

    def find_protein_context(self, hierarchy):
        return


class SiteList(object):
    def __init__(self, sites=None):
        """Class to hold information about multiple sites on a protein"""
        # List of Sites
        self.children = []
        self.parent = None
        self.id = None
        # Add Sites
        if sites: self.add_sites(sites)

    def add_sites(self, sites, update=True):
        if isinstance(sites, list):
            assert isinstance(sites[0], Site), 'Added sites must be of class Site'
            self.children.extend(sites)
        else:
            assert isinstance(sites, Site), 'Added sites must be of class Site'
            self.children.append(sites)
        if update:
            self.update()
        return self

    def update(self):
        # Number of events at the site
        self.num_sites = len(self.children)
        return self

    def summary(self):
        out = []
        out.append('SiteList Summary')
        out.append('{!s} Sites'.format(self.num_sites))
        for s in self.children:
            out.append('=====================================')
            out.append(s.summary())
        return '\n'.join(out)

    def sort(self, key, reverse=True):
        self.children.sort(key=key, reverse=True)
        return self

    def renumber(self, start_at=1):
        for i_c, c in enumerate(self.children): c.id = i_c + start_at
        return self


def events_table_to_events(event_table):
    # Only need the attributes: event.id, cluster.centroid, cluster.max

    events = []
    for record in event_table.itertuples():
        # print(record)

        # Spoof event id
        event_id = (record.dtag,
                    record.event_idx,
                    )

        # Spoof an empty cluster, then fix the centroid and max val
        event_cluster = PointCluster(id,
                                     points=[[0, 0, 0]],
                                     values=[0],
                                     )
        event_cluster.max = record.z_peak
        event_cluster.centroid = (record.x,
                                  record.y,
                                  record.z,
                                  )

        event = Event(id=event_id,
                      cluster=event_cluster,
                      )

        # centroids = flex.vec3_double([e.cluster.centroid
        #                               for e
        #                               in event.children
        #                               ]
        #                              )
        # event.info.centroid = centroids.mean()
        # event.info.approx_size = (centroids.min(), centroids.max())
        # event.info.num_events = len(event.children)
        #
        # print((event.id,
        #        event.cluster.centroid,
        #        event.cluster.max,
        #        )
        # )

        events.append(event)

    return events
