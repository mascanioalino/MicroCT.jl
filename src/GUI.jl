using Dash, DashCoreComponents, DashHtmlComponents, DashTable
export start_app
function start_app()
  app = dash()

  app.layout = html_div() do
      html_h1("Lung Dashboard"),
      html_div([
        html_div([
          html_h3("Enter directory path"),
          html_p("-If on Mac: go to directory, right click on the folder and hit option, then hit \"copy (name of file) as Pathname\""),
          html_p("-If on Windows: select the directory without going into it, hit \"Copy path\" on the tool bar on top"),
          html_p("Paste the path below")]),
        dcc_input(id = "path", type = "text", value = ""),
        html_div(id = "file_entered"),
        ]
      ),
      html_div([
        html_h3("Select function"),
        html_div([
          html_button("Create data csv", id="data"),
          ]),
        html_div(id = "data_creating"),
        html_div(id = "data_created"),
        html_div([
          html_button("Create histogram image and csv", id="hist", value="true")
          ]),
        html_div(id = "histogram_creating"),
        html_div(id = "histogram_created"),

        ]),
        # html_img(id = "histogram", src="/Users/mariaascanioalino/Dropbox (MIT)/microCT lung data/LaraTrainingData/Pre/AM27_Rec/histogram.png"),
      html_div([
        html_h3("Select option to do region of interest"),
        dcc_radioitems(
        options=[Dict("label"=> "Only lungs", "value"=> "roi"),
        Dict("label"=> "Whole mouse", "value"=> "no_roi")],
                        value="roi",
                        id="radio-options",),]
        ),
      html_div([
        html_h3("Select whether or not to save images"),
        dcc_checklist(
        options = [Dict("label"=> "Keep images", "value"=>"save_im")],
        id="save-images"
        )
      ])
  end

  #Entering file
  callback!(app, Output("file_entered", "children"), Input("path", "value")) do file
    html_div([html_p("You've entered $(file)")])
  end

  #Data collection
  callback!(app, Output("data_creating", "children"), Input("data", "n_clicks"), State("path", "value")) do n_clicks, file
    if length(file)!=0
      if isdir(file)
        html_div([html_p("Calculating data...")])

      else
        html_div([html_p("Please enter the correct path")])
      end
    else
      html_div([html_p("Please enter the correct path")])
    end
  end

  #Data collection
  callback!(app, Output("data_created", "children"), Input("data", "n_clicks"), State("path", "value")) do n_clicks, file
    if length(file)!=0
      if isdir(file)
        calculate_data(file)
        html_div([html_p("Done")])
      else
        html_div([html_p("")])
      end
    else
      html_div([html_p("")])
    end
  end

  #Show status
  callback!(app, Output("histogram_creating", "children"), Input("hist", "n_clicks"), State("path", "value")) do n_clicks, file
    if length(file)!=0
      if occursin(r"_Rec", file)
        html_div([html_p("Calculating histogram...")])
      else
        html_div([html_p("Calculating histograms...")])
      end
    else
      html_div([html_p("Please enter the correct path")])
    end
  end

  #Histogram creation
  callback!(app, Output("histogram_created", "children"), Input("hist", "n_clicks"), State("path", "value"), State("radio-options", "value"), State("save-images", "value")) do n_clicks, file, roi, save_im
    if length(file)!=0
      if occursin(r"_Rec", file)
        mouse = ReadStack(file)
        hist_data = calculate_singular_hist(mouse,roi=="roi",save_im=="save_im",file) # if it doesnt already exist
        html_div([html_p("Done")])
      elseif isdir(file)
        calculate_histograms(file,roi=="roi",save_im=="save_im")
        html_div([html_p("Done")])
      end
    else
      html_div([html_p("")])
    end
  end




  run_server(app, "0.0.0.0", 8080)
end
