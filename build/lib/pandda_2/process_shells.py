from collections import OrderedDict


class ProcessShells:
    def __init__(self,
                 process_shell,
                 processor,
                 ):
        self.process_shell = process_shell
        self.processor = processor

    def __call__(self,
                 shells,
                 reference,
                 grid,
                 tree,
                 ):
        args_list = []
        for shell_num, shell_dataset in shells.items():
            args = [shell_dataset,
                    reference,
                    grid,
                    tree,
                    shell_num,
                    ]
            args_list.append(args)

        self.processor(self.process_shell,
                       args_list,
                       )

    def __repr__(self):
        repr = OrderedDict()
        repr["process_shell"] = self.process_shell.repr()
        repr["processor"] = self.processor.repr()
        return repr
