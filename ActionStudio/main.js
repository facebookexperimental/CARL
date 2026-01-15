const express = require('express');

const app = express();
const port = 3000;

app.use(express.static('GithubPagesSite/docs'));

app.listen(port, () => {
    console.log(`Action Studio main server running at http://localhost:${port}`);
});
