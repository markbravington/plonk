## Turn OFF the goddamn byte-compiler
local({
	turn_off_byte_compiler <- function( ...) {
			# cat( 'TURNING OFF THE BLOODY BYTE-COMPILER, OR TRYING TO...\n')
			compiler:::enableJIT( 0)
			compiler:::compilePKGS( FALSE)
			# scatn( 'Success? %i', compiler:::enableJIT(-1))
		}
	environment( turn_off_byte_compiler) <- .GlobalEnv # avoid keeping exec-env of .First

	if( 'compiler' %in% loadedNamespaces()) {
		# cat( 'BYTE-COMPILER WAS ALREADY LOADED\n')
		turn_off_byte_compiler()
	} else {
		# did have 'attach' but surreptious sneaky auto-load happens...
		# ... possibly thru my own mistake! (now seemingly fixed)
		setHook( packageEvent( 'compiler', 'onLoad'), action='append',
				value=turn_off_byte_compiler)
	}
	# I would have thought the setHook above would kill the compiler, but it doesn't seem to do so properly...
	Sys.setenv( R_COMPILE_PKGS = '0')
	Sys.setenv( R_ENABLE_JIT = '0')
})
