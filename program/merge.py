import os
import shutil
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    # IO
    parser.add_argument("-o", "--old_pandda_path",
                        type=str,
                        help="The directory of the old pandda",
                        required=True
                        )
    parser.add_argument("-n", "--new_pandda_path",
                        type=str,
                        help="The directory of the new pandda (will be written to)",
                        required=True
                        )
    parser.add_argument("-m", "--merged_pandda_path",
                        type=str,
                        help="The directory of the new pandda (will be written to)",
                        required=True
                        )
    parser.add_argument("-ob", "--only_built",
                        type=bool,
                        help="Whether to only merge old events",
                        required=True
                        )

    args = parser.parse_args()

    return args


class Config:
    def __init__(self,
                 old_pandda_path,
                 new_pandda_path,
                 merged_pandda_path,
                 only_built,
                 ):
        self.old_pandda_path = old_pandda_path
        self.new_pandda_path = new_pandda_path
        self.merged_pandda_path = merged_pandda_path
        self.only_built = only_built


def get_config(args):
    config = Config(old_pandda_path=Path(args.old_pandda_path),
                    new_pandda_path=Path(args.new_pandda_path),
                    merged_pandda_path=Path(args.merged_pandda_path),
                    only_built=args.only_built,
                    )

    return config


class Output:
    def __init__(self, out_dir_path):
        self.merged_pandda_path = out_dir_path

    def attempt_mkdir(self, path):
        try:
            os.mkdir(str(path))
        except Exception as e:
            print(e)

    def attempt_remove(self, path):
        try:
            shutil.rmtree(path,
                          ignore_errors=True,
                          )
        except Exception as e:
            print(e)

    def make(self, overwrite=False):
        # Overwrite old results as appropriate
        if overwrite is True:
            self.attempt_remove(self.merged_pandda_path)

        # Make output dirs
        self.attempt_mkdir(self.merged_pandda_path)


def setup_output_directory(path, overwrite=False):
    output = Output(path)
    # output.make(overwrite)
    return output


class EventMapPath:
    def __init__(self):
        self.event_map_pattern = ""


class Event:
    def __init__(self,
                 row,
                 dtag,
                 event_idx,
                 occupancy,
                 analysed,
                 analysed_resolution,
                 cluster_size,
                 exclude_from_characterisation,
                 exclude_from_zmap_analysis,
                 global_correlation_to_average_map,
                 local_correlation_to_average_map,
                 map_uncertainty,
                 noisy_zmap,
                 r_free,
                 r_work,
                 refx,
                 refy,
                 refz,
                 rejected,
                 site_idx,
                 x,
                 y,
                 z,
                 z_mean,
                 z_peak,
                 interesting,
                 ligand_placed,
                 ligand_confidence,
                 comment,
                 viewed,
                 ):
        self.row = row
        self.dtag = dtag
        self.event_idx = event_idx
        self.occupancy = occupancy
        self.analysed = analysed
        self.analysed_resolution = analysed_resolution
        self.cluster_size = cluster_size
        self.exclude_from_characterisation = exclude_from_characterisation
        self.exclude_from_zmap_analysis = exclude_from_zmap_analysis
        self.global_correlation_to_average_map = global_correlation_to_average_map
        self.local_correlation_to_average_map = local_correlation_to_average_map
        self.map_uncertainty = map_uncertainty
        self.noisy_zmap = noisy_zmap
        self.r_free = r_free
        self.r_work = r_work
        self.refx = refx
        self.refy = refy
        self.refz = refz
        self.rejected = rejected
        self.site_idx = site_idx
        self.x = x
        self.y = y
        self.z = z
        self.z_mean = z_mean
        self.z_peak = z_peak
        self.interesting = interesting
        self.ligand_placed = ligand_placed
        self.ligand_confidence = ligand_confidence
        self.comment = comment
        self.viewed = viewed

    def to_record(self):
        record = {}

        record["dtag"] = self.dtag
        record["event_idx"] = self.event_idx
        record["1-BDC"] = self.occupancy
        record["analysed"] = self.analysed
        record["analysed_resolution"] = self.analysed_resolution
        record["cluster_size"] = self.cluster_size
        record["exclude_from_characterisation"] = self.exclude_from_characterisation
        record["exclude_from_zmap_analysis"] = self.exclude_from_zmap_analysis
        record["global_correlation_to_average_map"] = self.global_correlation_to_average_map
        record["local_correlation_to_average_map"] = self.local_correlation_to_average_map
        record["map_uncertainty"] = self.map_uncertainty
        record["noisy zmap"] = self.noisy_zmap
        record["r_free"] = self.r_free
        record["r_work"] = self.r_work
        record["refx"] = self.refx
        record["refy"] = self.refy
        record["refz"] = self.refz
        record["rejected - total"] = self.rejected
        record["site_idx"] = self.site_idx
        record["x"] = self.x
        record["y"] = self.y
        record["z"] = self.z
        record["z_mean"] = self.z_mean
        record["z_peak"] = self.z_peak
        record["Interesting"] = self.interesting
        record["Ligand Placed"] = self.ligand_placed
        record["Ligand Confidence"] = self.ligand_confidence
        record["Comment"] = self.comment
        record["Viewed"] = self.viewed

        return record

    @staticmethod
    def from_record(record):
        return Event(row=record,
                     dtag=record["dtag"],
                     event_idx=record["event_idx"],
                     occupancy=record["1-BDC"],
                     analysed=record["analysed"],
                     analysed_resolution=record["analysed_resolution"],
                     cluster_size=record["cluster_size"],
                     exclude_from_characterisation=record["exclude_from_characterisation"],
                     exclude_from_zmap_analysis=record["exclude_from_zmap_analysis"],
                     global_correlation_to_average_map=record["global_correlation_to_average_map"],
                     local_correlation_to_average_map=record["local_correlation_to_average_map"],
                     map_uncertainty=record["map_uncertainty"],
                     noisy_zmap=record["noisy zmap"],
                     r_free=record["r_free"],
                     r_work=record["r_work"],
                     refx=record["refx"],
                     refy=record["refy"],
                     refz=record["refz"],
                     rejected=record["rejected - total"],
                     site_idx=record["site_idx"],
                     x=record["x"],
                     y=record["y"],
                     z=record["z"],
                     z_mean=record["z_mean"],
                     z_peak=record["z_peak"],
                     interesting=record["Interesting"],
                     ligand_placed=record["Ligand Placed"],
                     ligand_confidence=record["Ligand Confidence"],
                     comment=record["Comment"],
                     viewed=record["Viewed"],
                     )

    def get_coords(self):
        return np.array([self.x, self.y, self.z])


class ModelPath:
    def __init__(self,
                 path,
                 ):
        self.path = path


class Site:
    def __init__(self,
                 site_idx,
                 events,
                 ):
        self.site_idx = site_idx
        event_coords = np.array([event.get_coords() for event in events])

        self.centroid = np.mean(event_coords,
                                axis=0,
                                )


def get_pandda_events(pandda_path):
    pandda_table_path = pandda_path / "analyses" / "pandda_inspect_events.csv"
    pandda_table = pd.read_csv(str(pandda_table_path))

    events = []
    for idx, row in pandda_table.iterrows():
        event = Event.from_record(row)
        events.append(event)

    return events


def get_closest_event(event,
                      old_pandda_events,
                      ):
    events_for_dtag = []
    for old_event in old_pandda_events:
        if old_event.dtag == event.dtag:
            events_for_dtag.append(old_event)

    if len(events_for_dtag) == 0:
        return None, 0
    else:
        distances = []
        for dtag_event in events_for_dtag:
            distance = get_distance(event.get_coords(),
                                    dtag_event.get_coords(),
                                    )
            distances.append(distance)

        closest_event_idx = np.argmin(distances)

        closest_event = events_for_dtag[closest_event_idx]
        distance = distances[closest_event_idx]

        return closest_event, distance


def get_events_to_merge(old_pandda_events,
                        new_pandda_events,
                        only_built,
                        distance_cutoff=6.0,
                        ):
    matches_old_to_new = {}
    unmatched_events = []
    # Get closest old event to the new events
    for event in new_pandda_events:
        closest_event, distance = get_closest_event(event,
                                                    old_pandda_events,
                                                    )

        if closest_event is None:
            unmatched_events.append(event)
        elif distance > distance_cutoff:
            unmatched_events.append(event)
        else:
            matches_old_to_new[(closest_event.dtag, closest_event.event_idx)] = event

    return unmatched_events, matches_old_to_new


def copy_file(path_1,
              path_2,
              ):
    shutil.copy(str(path_1),
                str(path_2),
                )


def copy_sites_table(old_pandda_path,
                     new_pandda_path,
                     ):
    old_pandda_site_table_path = old_pandda_path / "analyses" / "pandda_inspect_sites.csv"
    new_pandda_site_table_path = new_pandda_path / "analyses" / "pandda_inspect_sites.csv"
    copy_file(old_pandda_site_table_path,
              new_pandda_site_table_path,
              )


def get_sites(old_pandda_events):
    event_sites = {}

    for event in old_pandda_events:
        if event.site_idx in event_sites:
            event_sites[event.site_idx].append(event)
        else:
            event_sites[event.site_idx] = [event]

    sites = {}
    for site_idx in event_sites:
        site = Site(site_idx,
                    event_sites[site_idx],
                    )
        sites[site_idx] = site

    return sites


def get_distance(point_1, point_2):
    return np.linalg.norm(point_1 - point_2)


def get_closest_site(event,
                     sites,
                     ):
    site_idxs = []
    distances = []
    for site_idx, site in sites.items():
        distance = get_distance(event.get_coords(),
                                site.centroid,
                                )
        site_idxs.append(site_idx)
        distances.append(distance)

    closest_site_idx = site_idxs[np.argmin(distances)]
    return sites[closest_site_idx]


def update_events(unmatched_events,
                  matches_old_to_new,
                  old_pandda_events,
                  sites,
                  ):
    merged_events = {}
    event_mapping = {}

    # Get base events
    for event in old_pandda_events:
        merged_events[(event.dtag, event.event_idx)] = event
        # event_mapping[(event.dtag, event.event_idx)] = event

    # Update for merged events
    for old_event_idx, new_event in matches_old_to_new.items():
        old_event = merged_events[old_event_idx]
        closest_site = merged_events[old_event_idx].site_idx
        new_event = Event(row=new_event.row,
                          dtag=old_event.dtag,
                          event_idx=new_event.event_idx,
                          occupancy=new_event.occupancy,
                          analysed=new_event.analysed,
                          analysed_resolution=new_event.analysed_resolution,
                          cluster_size=new_event.cluster_size,
                          exclude_from_characterisation=new_event.exclude_from_characterisation,
                          exclude_from_zmap_analysis=new_event.exclude_from_zmap_analysis,
                          global_correlation_to_average_map=new_event.global_correlation_to_average_map,
                          local_correlation_to_average_map=new_event.local_correlation_to_average_map,
                          map_uncertainty=new_event.map_uncertainty,
                          noisy_zmap=new_event.noisy_zmap,
                          r_free=new_event.r_free,
                          r_work=new_event.r_work,
                          refx=new_event.refx,
                          refy=new_event.refy,
                          refz=new_event.refz,
                          rejected=new_event.rejected,
                          site_idx=closest_site,
                          x=new_event.x,
                          y=new_event.y,
                          z=new_event.z,
                          z_mean=new_event.z_mean,
                          z_peak=new_event.z_peak,
                          interesting=new_event.interesting,
                          ligand_placed=new_event.ligand_placed,
                          ligand_confidence=new_event.ligand_confidence,
                          comment=new_event.comment,
                          viewed=new_event.viewed,
                          )
        merged_events[old_event_idx] = new_event
        event_mapping[old_event_idx] = new_event

    # Update for completely new events
    for unmatched_event in unmatched_events:
        dataset_events = [event_idx for event_idx in merged_events.keys() if event_idx[0] == unmatched_event.dtag]
        if len(dataset_events) > 0:
            new_event_idx = len(dataset_events) + 1
        else:
            new_event_idx = 1

        closest_site = get_closest_site(unmatched_event,
                                        sites,
                                        )

        new_event = Event(row=unmatched_event.row,
                          dtag=unmatched_event.dtag,
                          event_idx=new_event_idx,
                          occupancy=unmatched_event.occupancy,
                          analysed=unmatched_event.analysed,
                          analysed_resolution=unmatched_event.analysed_resolution,
                          cluster_size=unmatched_event.cluster_size,
                          exclude_from_characterisation=unmatched_event.exclude_from_characterisation,
                          exclude_from_zmap_analysis=unmatched_event.exclude_from_zmap_analysis,
                          global_correlation_to_average_map=unmatched_event.global_correlation_to_average_map,
                          local_correlation_to_average_map=unmatched_event.local_correlation_to_average_map,
                          map_uncertainty=unmatched_event.map_uncertainty,
                          noisy_zmap=unmatched_event.noisy_zmap,
                          r_free=unmatched_event.r_free,
                          r_work=unmatched_event.r_work,
                          refx=unmatched_event.refx,
                          refy=unmatched_event.refy,
                          refz=unmatched_event.refz,
                          rejected=unmatched_event.rejected,
                          site_idx=closest_site.site_idx,
                          x=unmatched_event.x,
                          y=unmatched_event.y,
                          z=unmatched_event.z,
                          z_mean=unmatched_event.z_mean,
                          z_peak=unmatched_event.z_peak,
                          interesting=unmatched_event.interesting,
                          ligand_placed=unmatched_event.ligand_placed,
                          ligand_confidence=unmatched_event.ligand_confidence,
                          comment=unmatched_event.comment,
                          viewed=unmatched_event.viewed,
                          )

        merged_events[(unmatched_event.dtag, new_event_idx)] = new_event
        event_mapping[(unmatched_event.dtag, new_event_idx)] = unmatched_event

    return merged_events, event_mapping

    #         # Get event dict
    # event_dict = {}
    # for event in old_pandda_events:
    #     if event.dtag not in event_dict:
    #         event_dict[event.dtag] = {}
    #
    #     event_dict[event.dtag][event.event_idx] = event
    #
    # # Create new events
    # event_mapping = {}
    # new_events = []
    # for event in events_to_merge:
    #     if event.dtag in event_dict:
    #         new_event_idx = max(event_dict[event.dtag].keys()) + 1
    #     else:
    #         event_dict[event.dtag] = {}
    #         new_event_idx = 1
    #
    #     closest_site = get_closest_site(event,
    #                                     sites,
    #                                     )
    #     new_event = Event(row=event.row,
    #                       dtag=event.dtag,
    #                       event_idx=new_event_idx,
    #                       occupancy=event.occupancy,
    #                       analysed=event.analysed,
    #                       analysed_resolution=event.analysed_resolution,
    #                       cluster_size=event.cluster_size,
    #                       exclude_from_characterisation=event.exclude_from_characterisation,
    #                       exclude_from_zmap_analysis=event.exclude_from_zmap_analysis,
    #                       global_correlation_to_average_map=event.global_correlation_to_average_map,
    #                       local_correlation_to_average_map=event.local_correlation_to_average_map,
    #                       map_uncertainty=event.map_uncertainty,
    #                       noisy_zmap=event.noisy_zmap,
    #                       r_free=event.r_free,
    #                       r_work=event.r_work,
    #                       refx=event.refx,
    #                       refy=event.refy,
    #                       refz=event.refz,
    #                       rejected=event.rejected,
    #                       site_idx=closest_site.site_idx,
    #                       x=event.x,
    #                       y=event.y,
    #                       z=event.z,
    #                       z_mean=event.z_mean,
    #                       z_peak=event.z_peak,
    #                       interesting=event.interesting,
    #                       ligand_placed=event.ligand_placed,
    #                       ligand_confidence=event.ligand_confidence,
    #                       comment=event.comment,
    #                       viewed=event.viewed,
    #                       )
    #     new_events.append(new_event)
    #     # Map FINAL event id to ORIGINAL event id
    #     event_mapping[(new_event.dtag, new_event.event_idx)] = event
    #
    # final_events = new_events + old_pandda_events

    # return merged_events, event_mapping


def copy_models(events_to_merge,
                old_pandda_path,
                new_pandda_path,
                ):
    for event in events_to_merge:
        old_model_path = old_pandda_path / "processed_datasets" / event.dtag / "modelled_structures" / "{}-pandda-model.pdb".format(
            event.dtag)
        new_model_path = new_pandda_path / "processed_datasets" / event.dtag / "modelled_structures" / "{}-pandda-model.pdb".format(
            event.dtag)
        copy_file(old_model_path,
                  new_model_path,
                  )


def copy_directory(directory_1, directory_2):
    shutil.copytree(str(directory_1),
                    str(directory_2),
                    )


def sync_event_dirs(final_events,
                    new_pandda_path,
                    merged_pandda_path,
                    event_mapping_merged_to_original,
                    ):
    pandda_event_map_pattern = "{dtag}-event_{event_idx}_1-BDC_{occupancy}_map.ccp4"
    pandda_event_model_pattern = "{dtag}-pandda-model.pdb"

    sync_records = {}
    for event_id, event in final_events.items():
        # If event is not to be updated, skip
        if event_id not in event_mapping_merged_to_original:
            sync_records[(event.dtag, event.event_idx)] = "\t\tNo changes to event {}".format(event_id)
            continue

        # If dataset is completely new, copy from new pandda and continue
        merged_event_dir_path = merged_pandda_path / "processed_datasets" / event.dtag
        new_event_dir_path = new_pandda_path / "processed_datasets" / event.dtag
        if not merged_event_dir_path.exists():
            copy_directory(new_event_dir_path,
                           merged_event_dir_path,
                           )
            sync_records[(event.dtag, event.event_idx)] = "New dataset only analysed in new pandda"
            continue

        # If the event is new or updated, copy the new event map and model if they are there
        new_event = event_mapping_merged_to_original[event_id]

        # Copy event map
        new_event_map_path = new_event_dir_path / pandda_event_map_pattern.format(dtag=new_event.dtag,
                                                                                  event_idx=new_event.event_idx,
                                                                                  occupancy=new_event.occupancy,
                                                                                  )
        merged_event_map_path = merged_event_dir_path / pandda_event_map_pattern.format(dtag=event.dtag,
                                                                                        event_idx=event.event_idx,
                                                                                        occupancy=event.occupancy,
                                                                                        )
        if new_event_map_path.exists():
            copy_file(new_event_map_path,
                      merged_event_map_path,
                      )
            sync_records[(event.dtag, event.event_idx)] = "Copied event map for new event"
        else:
            sync_records[(event.dtag, event.event_idx)] = "No event map for this event"

        # Copy model
        new_model_path = new_event_dir_path / "modelled_structures" / pandda_event_model_pattern.format(
            dtag=new_event.dtag,
            event_idx=new_event.event_idx,
            occupancy=new_event.occupancy,
        )
        merged_model_path = merged_event_dir_path / "modelled_structures" / pandda_event_model_pattern.format(
            dtag=event.dtag,
            event_idx=event.event_idx,
            occupancy=event.occupancy,
        )

        if new_model_path.exists():
            copy_file(new_model_path,
                      merged_model_path,
                      )
            sync_records[(event.dtag, event.event_idx)] = "Copied model for new event"
        else:
            sync_records[(event.dtag, event.event_idx)] = "No model for this event"

    return sync_records


def make_event_table(final_events):
    records = []
    for event_idx, event in final_events.items():
        record = event.to_record()
        records.append(record)

    event_table_merged = pd.DataFrame(records)

    return event_table_merged


def output_events_table(event_table,
                        pandda_path,
                        ):
    event_table_analyse_path = pandda_path / "analyses" / "pandda_inspect_events.csv"
    event_table.to_csv(str(event_table_analyse_path))
    event_table_inspect_path = pandda_path / "analyses" / "pandda_analyse_events.csv"
    event_table.to_csv(str(event_table_inspect_path))


def main():
    args = parse_args()

    config = get_config(args)

    output = setup_output_directory(config.merged_pandda_path)

    print("Copying old pandda directory...")
    if not output.merged_pandda_path.exists():
        copy_directory(config.old_pandda_path,
                       output.merged_pandda_path,
                       )
    print("\tCoptied old pandda directory!")

    print("Getting old pandda events...")
    old_pandda_events = get_pandda_events(config.old_pandda_path)
    print("\tGot {} events from the old pandda".format(len(old_pandda_events)))

    print("Getting new pandda events...")
    new_pandda_events = get_pandda_events(config.new_pandda_path)
    print("\tGot {} events from the new pandda".format(len(new_pandda_events)))

    print("Figuring out which events to merge...")
    unmatched_events, matches_old_to_new = get_events_to_merge(old_pandda_events,
                                                               new_pandda_events,
                                                               config.only_built,
                                                               )
    print("\tGot {} unmatched events...".format(len(unmatched_events)))
    print("\tGot {} matched events...".format(len(matches_old_to_new)))

    print("Getting sites of old pandda")
    sites = get_sites(old_pandda_events)
    print("\tGot {} sites".format(len(sites)))

    print("Updating events...")
    final_events, event_mapping = update_events(unmatched_events,
                                                matches_old_to_new,
                                                old_pandda_events,
                                                sites,
                                                )
    print("\tNow have {} final events".format(len(final_events)))

    # Update the names of event maps
    print("Syncing new systems and their event maps...")
    sync_records = sync_event_dirs(final_events,
                                   config.new_pandda_path,
                                   output.merged_pandda_path,
                                   event_mapping,
                                   )
    print("\tSynced event dirs!")
    for key, record in sync_records.items():
        print("\t{}: {}".format(key, record))

    print("Making event table of final events...")
    events_table = make_event_table(final_events)
    print("\tMade event table!")

    print("Outputting event table to {}...".format(output.merged_pandda_path))
    output_events_table(events_table,
                        output.merged_pandda_path,
                        )
    print("\tOutput event table!")


if __name__ == "__main__":
    main()
