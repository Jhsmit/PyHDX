from dask.distributed import as_completed
import param

import panel as pn



class ASyncProgressBar(param.Parameterized):
    completed = param.Integer(default=0)
    num_tasks = param.Integer(default=10, bounds=(1, None))

    async def run(self, futures):
        async for task in as_completed(futures):
            with pn.io.unlocked():
                self.completed += 1

    @property
    def value(self):
        return int(100 * (self.completed / self.num_tasks))

    def reset(self):
        self.completed = 0

    def increment(self):
        self.completed += 1

    @param.depends('completed', 'num_tasks')
    def view(self):
        if self.value:
            return pn.widgets.Progress(active=True, value=self.value, align="center", sizing_mode="stretch_width")
        else:
            return None
