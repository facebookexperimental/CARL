import React from 'react';
import { createRoot } from 'react-dom/client';
import App from './App';

/*const Greet = () => {
    const canvasRef = React.useRef(null);
    const buttonRef = React.useRef(null);
    const examples = new Set();
    const definitions = new Set();
    const [examplesCount, setExamplesCount] = React.useState(0);
    const [definitionsCount, setDefinitionsCount] = React.useState(0);
    async function bootstrapAsync() {
        const carl = await CarlIntegration.CreateAsync();
        carl.onExampleCreated = example => {
            examples.add(example);
            setExamplesCount(examples.size);
            example.onDisposed = () => {
                examples.delete(example);
                setExamplesCount(examples.size);
            };
        };
        carl.onDefinitionCreated = definition => {
            definitions.add(definition);
            setDefinitionsCount(definitions.size);
            definition.onDisposed = () => {
                definitions.delete(definition);
                setDefinitionsCount(definitions.size);
            };
        };
        const immersiveExperience = await initializeImmersiveExperienceAsync(canvasRef.current, carl);
        buttonRef.current.addEventListener("click", () => {
            immersiveExperience.enterImmersiveMode();
        });
    }
    React.useEffect(() => {
        bootstrapAsync();
    }, []);

    return <>
        <h1>Hello, world!</h1>
        <p>How's life?</p>
        <p>You have {examplesCount} examples.</p>
        <p>You have {definitionsCount} definitions.</p>
        <button ref={buttonRef}>Enter XR</button>
        <canvas ref={canvasRef} width={16} height={16} style={{ display: 'none' }}></canvas>
    </>;
};*/

const container = document.getElementById('root');
const root = createRoot(container);
root.render(<App />);

