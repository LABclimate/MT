
import xarray as xr

def concat(output, input):
    try:    return(xr.concat([output, input], dim='time'))
    except: return(input)
    
