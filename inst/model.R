
      library( TMB )
      dyn.load( dynlib('inst/model') )
      setwd('C:/Users/test/R_projects/FishMap')
      load( 'All_inputs.RData' )
      Obj = MakeADFun(data=All_inputs[['data']], parameters=All_inputs[['parameters']], random=All_inputs[['random']], All_inputs[['Other_inputs']])
      save( Obj, file='Obj.RData')
    