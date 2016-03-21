
function sfileinfo, fn
    info = file_info(fn)
    case !version.os_family of
        'Windows': begin
            spawn, 'powershell $a = Get-Date -Date ((ls '+fn+$
                ').LastWriteTimeUtc) -UFormat %s;$a', mtime
            mtime = double(mtime[0]) & info.mtime = mtime & end
        else:
    endcase
    return, info
end

fn = shomedir()+'/thg_l1_ast_fsmi_20130501_v01.cdf'
fn = shomedir()+'/test.txt'
info = sfileinfo(fn)
print, time_string(info.mtime)
end