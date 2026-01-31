import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import { initializeImmersiveExperienceAsync } from 'carl-actionstudio-immersiveexperience';
import { initializeNativeIntegrationAsync } from 'carl-actionstudio-nativeintegration';

const Greet = () => {
    const canvasRef = React.useRef(null);
    async function bootstrapAsync() {
        const carl = await initializeNativeIntegrationAsync();
        const inputSample = carl.createInputSample();

        let ipr = null;
        await initializeImmersiveExperienceAsync(canvasRef.current, inputSample, sample => {
            if (ipr === null) {
                ipr = new carl.InProgressRecording();
                setTimeout(() => {
                    const recording = new carl.Recording(ipr);
                    ipr.delete();
                    //ipr = null;

                    const inspector = new carl.RecordingInspector(recording);
                    const startT = inspector.startTimestamp();
                    const endT = inspector.endTimestamp();
                    inspector.delete();

                    const example = new carl.Example(recording, startT, endT);
                    recording.delete();

                    const bytes = example.serialize();
                    example.delete();

                    const jsBytes = new Uint8Array(bytes.size());
                    for (let idx = 0; idx < jsBytes.length; ++idx) {
                        jsBytes[idx] = bytes.get(idx);
                    }
                    bytes.delete();

                    const buffer = jsBytes.buffer;
                    const blob = new Blob([buffer], { type: "application/octet-stream" });
                    const link = document.createElement('a');
                    link.href = URL.createObjectURL(blob);
                    link.download = "example_0.bin";
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    URL.revokeObjectURL(link.href);
                }, 5000);
            }
            if (!ipr.isDeleted()) {
                ipr.addInput(sample);
            }
        });
        console.log("Done!");
    }
    React.useEffect(() => {
        bootstrapAsync();
    }, []);

    return <>
        <h1>Hello, world!</h1>
        <p>How's life?</p>
        <canvas ref={canvasRef} width={400} height={400} style={{ border: '1px solid black' }}></canvas>
    </>;
};

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(<Greet />);
