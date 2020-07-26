import param
from panel.widgets.input import Widget, _BkTextInput, StaticText
from panel.widgets import Spinner


class NumericInput(Widget):
    value = param.Number(default=0, allow_None=True, bounds=[None, None])

    placeholder = param.Number(default=None)

    start = param.Number(default=None, allow_None=True)

    end = param.Number(default=None, allow_None=True)

    _widget_type = _BkTextInput

    formatter = param.Parameter(default=None)

    _rename = {'name': 'title', 'formatter': None, 'start': None, 'end': None}

    def __init__(self, **params):
        if params.get('value') is None:
            value = params.get('start', self.value)
            if value is not None:
                params['value'] = value
        super(NumericInput, self).__init__(**params)

    def _bound_value(self, value):
        if self.start is not None:
            value = max(value, self.start)
        if self.end is not None:
            value = min(value, self.end)
        return value

    def _format_value(self, value):
        if self.formatter is not None:
            value = self.formatter.format(value)
        else:
            value = str(value)
        return value

    def _process_param_change(self, msg):
        msg.pop('formatter', None)

        if 'start' in msg:
            start = msg.pop('start')
            self.param.value.bounds[0] = start
        if 'end' in msg:
            end = msg.pop('end')
            self.param.value.bounds[1] = end

        if 'value' in msg and msg['value'] is not None:
            msg['value'] = self._format_value(self.value)
        if 'placeholder' in msg and msg['placeholder'] is not None:
            msg['placeholder'] = self._format_value(self.placeholder)
        return msg

    def _process_property_change(self, msg):
        if 'value' in msg and msg['value'] is not None:
            try:
                value = float(msg['value'])
                msg['value'] = self._bound_value(value)
                if msg['value'] != value:
                    self.param.trigger('value')
            except ValueError:
                msg.pop('value')
        if 'placeholder' in msg and msg['placeholder'] is not None:
            try:
                msg['placeholder'] = self._format_value(float(msg['placeholder']))
            except ValueError:
                msg.pop('placeholder')
        return msg


class ColoredStaticText(StaticText):
    _format = '<b>{title}</b>: <span class ="panel-colored-statictext">{value}</span>'
