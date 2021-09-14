import Pkg; import FileWatching: watch_file


function filewatcher(filename)
    while true
        sleep(0.1)
        event = watch_file("./")
        if event.changed
            println("File changed, running")
            try
                include(filename)
            catch err
                showerror(stdout, err, catch_backtrace())
            end
        end
    end
end