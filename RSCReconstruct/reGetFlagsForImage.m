function flags=reGetFlagsForImage(emFlags,iImg)
switch emFlags.type
    case 'logical'
        flags=emFlags.flags(:,:,:,iImg);
    otherwise
        disp(['Unknown flag type: ' emFlags.type]);
end;
end

