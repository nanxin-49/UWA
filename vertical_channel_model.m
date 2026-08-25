function output = vertical_channel_model(paramsV)
%VERTICAL_CHANNEL_MODEL Public upward vertical acoustic channel API.
%   The implementation lives under src/. This root wrapper preserves the
%   established call signature while initializing the project path.

setup_vertical_project();
output = vertical_channel_model_impl(paramsV);
end
