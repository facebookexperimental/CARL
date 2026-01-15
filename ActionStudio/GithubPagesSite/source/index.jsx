import * as React from 'react';
import * as ReactDOM from 'react-dom/client';
import { initializeImmersiveExperience } from 'carl-actionstudio-immersiveexperience';
import { testNativeIntegration } from 'carl-actionstudio-nativeintegration';

const Greet = () => {
    const canvasRef = React.useRef(null);
    React.useEffect(() => {
        initializeImmersiveExperience(canvasRef.current);
        testNativeIntegration();
    }, []);

    return <>
        <h1>Hello, world!</h1>
        <p>How's life?</p>
        <canvas ref={canvasRef} width={400} height={400} style={{ border: '1px solid black' }}></canvas>
    </>;
};

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(<Greet />);
