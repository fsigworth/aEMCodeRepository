function flags=rePutFlagsForImage(flags,iImg)
switch emFlags.type
    case 'logical'
        emFlags.flags(:,:,:,iImg)=flags;
    otherwise
        disp(['Unknown flag type: ' emFlags.type]);
end;
end

