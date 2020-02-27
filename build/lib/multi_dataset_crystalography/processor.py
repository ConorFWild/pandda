class ProcessModelSeriel:

    def __init__(self):
        pass

    def __call__(self, model):
        model.process()

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

class ProcessModelQsub:
    def __init__(self):
        self.script_paths = []
        self.script_names = []

    def __call__(self, model):

        # Define submission script
        submission_script = ""

        # Record script name

    def __enter__(self):
        # Define python script

        import_pickle = "" # Import pickle

        load_model = ""# Load model

        call_model = ""# call model

        # Save script to file

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Wait for script to complete

        # Resubmit any scripts that were not run

        # Delete python script

        # Delete runscripts

        pass