class PanDDARunLog:
    # Produces a log of data processing from dataset, model, graphs and HTML

    def __init__(log):
        log(pandda_dataset)

        for model in pandda_models: log(model)

        log(pandda_run_graphs)

        log(pandda_run_html)