export magmaInit, magmaFinalize

# magma_init
function magmaInit()
	ccall((:magma_init, libmagma),Cint,())
end

# magma_finalize
function magmaFinalize()
	ccall((:magma_finalize, libmagma),Cint,())
end
