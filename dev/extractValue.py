import panel as pn
import param
from panel.pane import HTML, Markdown

class testclass(HTML):

    def __init__(self,**params):
        super().__init__(**params)
        self.object =\
            """ 
                <!DOCTYPE html>
                <html>
                <body>
                
                <h2>JavaScript Class</h2>
                
                <p>How to use a JavaScript Class.</p>
                
                <p id="demo"></p>
                
                <script>
                class Car {
                  constructor(name, year) {
                    this.name = name;
                    this.year = year;
                  }
                }
                
                myCar = new Car("Ford", 2014);
                document.getElementById("demo").innerHTML =
                myCar.name + " " + myCar.year;
                </script>
                
                </body>
                </html>
        """
test = testclass(sizing_mode='stretch_both')
test.servable()
